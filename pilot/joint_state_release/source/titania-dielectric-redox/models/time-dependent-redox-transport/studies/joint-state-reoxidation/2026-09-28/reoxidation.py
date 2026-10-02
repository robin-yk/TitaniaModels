"""R600 measured reoxidation demand applied to saved joint-state reduction profiles.

No transport/energy fitting. One measured CO removes one vacancy and two electrons.
The finite-volume residual and microscopic free energy are imported unchanged.
"""

from pathlib import Path as _TitaniaPath
_TITANIA_PROJECT_ROOT = next(p for p in _TitaniaPath(__file__).resolve().parents if p.name == 'titania-dielectric-redox')
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
from types import ModuleType

import numpy as np
from scipy.optimize import brentq
from scipy.linalg import LinAlgWarning

HERE = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28'
PROJECT = _TITANIA_PROJECT_ROOT
PREVIOUS = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28'
OLD = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/versions/0.1.0'
SOURCE = _TITANIA_PROJECT_ROOT / 'models/atomic-layer-and-bulk-defect-distribution/studies/site-population/reoxidation-source.csv'
sys.path.insert(0, str(PREVIOUS))
import connection as c


def measured():
    data = np.genfromtxt(SOURCE, delimiter=',', names=True)
    te = np.r_[0., data['time_s'], 600.]
    re = np.r_[data['rate_umol_g_s'][0], data['rate_umol_g_s'],
               data['rate_umol_g_s'][-1]]
    raw = np.r_[0., np.cumsum(np.diff(te) * (re[1:] + re[:-1]) / 2)]
    cumulative = np.maximum.accumulate(raw)
    # Piecewise-linear cumulative demand: each interval uses its measured mean.
    def amount(t):
        return float(np.interp(t, te, cumulative))
    def rate(t):
        j = min(np.searchsorted(te, t, side='right') - 1, len(te) - 2)
        return float((cumulative[j+1] - cumulative[j]) / (te[j+1] - te[j]))
    return amount, rate, dict(time_s=te.tolist(), measured_rate_umol_g_s=re.tolist(),
                             raw_cumulative_umol_g=raw.tolist(),
                             imposed_cumulative_umol_g=cumulative.tolist(),
                             processing='Constant endpoints to 0 and 600 s; trapezoidal integration; '
                             'running maximum of cumulative CO; linear cumulative interpolation; no rescaling.')


def engine(delta, binding):
    adapter = ModuleType('pilot')
    adapter.HERE, adapter.T, adapter.pt = OLD, c.T, c.pt
    adapter.build = lambda n: c.build(n, delta, binding)
    adapter.quantities, adapter.source = c.quantities, measured
    sys.modules['pilot'] = adapter
    spec = importlib.util.spec_from_file_location('original_transport', _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/versions/0.1.0/transient.py')
    e = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(e)
    return e


def restore(tr, checkpoint):
    """Recover the omitted electron multiplier without equilibrating vacancies."""
    L = tr.L
    Y = np.r_[np.array(checkpoint['vacancy_chemical_potential_eV']) / tr.kT,
              0., np.array(checkpoint['phi_eV_per_e']) / tr.kT]
    target = checkpoint['electron_total_umol_g']
    def error(ye):
        Y[L] = ye
        return tr.state(Y, False)[0][L] - target
    Y[L] = brentq(error, -1000., 1000., xtol=1e-13)
    g, x, q = tr.state(Y, False)
    row = c.quantities(tr.m, tr.lay, x, Y[L+1:] * tr.kT)
    checks = {}
    for name in ('vacancy_umol_g', 'electron_umol_g'):
        checks['restart_' + name + '_max_error'] = float(np.max(np.abs(
            np.array(row[name]) - np.array(checkpoint[name]))))
        assert checks['restart_' + name + '_max_error'] < 1e-7
    checks['restart_poisson_residual'] = float(np.max(np.abs(g[L+1:])))
    return Y, g[:L].copy(), x, q, checks


def run(args):
    initial_path = PREVIOUS / args.initial
    previous = json.loads(initial_path.read_text())
    p = previous['assumed_parameters']
    e = engine(p['delta_eV'], p['binding_eV'])
    tr = e.Transport(previous['n_bulk'], previous['D_nm2_s'])
    tr.null_tail = True
    saved = previous['checkpoints'][-1]
    assert saved['time_s'] == 1800.
    Y, old, x, q, checks = restore(tr, saved)
    checks['initial_jacobian_relative_error'] = e.jacobian_check(tr, Y)
    amount, rate, source = measured()
    N0 = float(old.sum())
    neutral0 = float(Y[0] + 2 * Y[tr.L])
    checks.update(max_local_balance_umol_g=0., max_vacancy_balance_umol_g=0.,
                  max_electron_balance_umol_g=0., min_transport_dissipation=1e100,
                  max_open_free_energy_inequality_umol_eV_g=-1e100)
    rows, log, retries = [], [], []
    t, dt = 0., .001
    Fprevious = tr.free_energy(x, q)
    start = time.monotonic()
    last_dt, last_it, last_err = 0., 0, 0.

    def row():
        r = tr.row(Y, x, q, t, N0 - amount(t), last_dt, last_it, last_err)
        r.update(withdrawn_CO_umol_g=amount(t), demand_umol_g_s=rate(t),
                 neutral_surface_mu_change_kT=float(Y[0] + 2*Y[tr.L] - neutral0),
                 full_state_Y=Y.tolist(), outward_first_face_flux_umol_g_s=-r['inward_flux_umol_g_s'][0])
        return r
    rows.append(row())
    stops = sorted(set([.01, .1, .5, 1., 2., 5., 10., 15., 20., 25., 30., 40.,
                        60., 90., 120., 180., 300., 450., 600., args.end] +
                       source['time_s']))
    status = 'completed'
    failure = None
    while t < args.end - 1e-9:
        stop = min(s for s in stops if s > t + 1e-9)
        dt = min(dt, args.max_step, stop-t, args.end-t)
        target = N0 - amount(t+dt)
        added = amount(t) - amount(t+dt)
        try:
            Yn, xn, qn, J, it, err = tr.step(Y.copy(), old, added, dt, target)
        except (RuntimeError, LinAlgWarning) as exc:
            retries.append(dict(t_s=t, dt_s=dt, error=str(exc)))
            dt *= .5
            if dt < args.min_step:
                status = 'numerical continuation stopped; physical cause requires diagnosis'
                failure = dict(time_s=t, attempted_dt_s=2*dt, error=str(exc))
                break
            continue
        g = tr.state(Yn, False)[0]
        local = g[:tr.L] - old
        local[0] -= added
        local[:-1] += dt*J
        local[1:] -= dt*J
        diss = float(tr.kT * J @ (Yn[:tr.L-1] - Yn[1:tr.L]))
        F = tr.free_energy(xn, qn)
        # Convex free energy: ΔF <= ΔN (μV,s + 2 μe) - dt*dissipation.
        inequality = F - Fprevious - added * tr.kT * (Yn[0]+2*Yn[tr.L]) + dt*diss
        checks['max_local_balance_umol_g'] = max(checks['max_local_balance_umol_g'], float(np.max(np.abs(local))))
        checks['max_vacancy_balance_umol_g'] = max(checks['max_vacancy_balance_umol_g'], abs(float(g[:tr.L].sum()-target)))
        checks['max_electron_balance_umol_g'] = max(checks['max_electron_balance_umol_g'], abs(float(g[tr.L]-2*target)))
        checks['min_transport_dissipation'] = min(checks['min_transport_dissipation'], diss)
        checks['max_open_free_energy_inequality_umol_eV_g'] = max(checks['max_open_free_energy_inequality_umol_eV_g'], inequality)
        Y, old, x, q, t, Fprevious = Yn, g[:tr.L].copy(), xn, qn, t+dt, F
        last_dt, last_it, last_err = dt, it, err
        log.append(dict(t_s=t, dt_s=dt, iterations=it, residual=err,
                        theta_BRI_pct=100*g[0]/tr.lay.c_bri,
                        neutral_mu_change_kT=float(Y[0]+2*Y[tr.L]-neutral0),
                        outward_flux_umol_g_s=float(-J[0])))
        if abs(t-stop) < 1e-8:
            r = row()
            rows.append(r)
            if t in (.1, 1., 5., 10., 15., 20., 30., 60., 120., 300., 600.):
                print(json.dumps(dict(t=t, theta=r['theta_BRI_pct'], CO=amount(t),
                      ds=r['neutral_surface_mu_change_kT'], elapsed=time.monotonic()-start)), flush=True)
        dt = min(args.max_step, dt*1.7)
    if rows[-1]['time_s'] != t:
        rows.append(row())
    checks['final_jacobian_relative_error'] = e.jacobian_check(tr, Y)
    assert checks['max_vacancy_balance_umol_g'] < 1e-7
    assert checks['max_electron_balance_umol_g'] < 1e-7
    assert checks['min_transport_dissipation'] > -1e-8
    assert checks['max_open_free_energy_inequality_umol_eV_g'] < 1e-6
    files = [Path(__file__), initial_path, _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/connection.py', _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/versions/0.1.0/transient.py', SOURCE]
    output = dict(model='model-tio2-surface-bulk-joint-state v0.1.0, measured-demand reoxidation',
                  starting_profile=str(initial_path.relative_to(PROJECT)),
                  parameters=p, T_K=c.T, diameter_nm=c.DIAMETER, n_bulk=tr.L-12,
                  D_nm2_s=tr.D, max_step_s=args.max_step, min_step_s=args.min_step,
                  start_inventory_umol_g=N0, end_time_s=t, status=status, failure=failure,
                  source=source, checks=checks, checkpoints=rows, steps=log, retries=retries,
                  elapsed_s=time.monotonic()-start,
                  sha256={str(f.relative_to(PROJECT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    path = HERE / args.output
    path.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(dict(output=str(path), status=status, end_s=t, checks=checks,
                          accepted=len(log), rejected=len(retries))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--initial', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-step', type=float, default=2.)
    parser.add_argument('--min-step', type=float, default=1e-7)
    parser.add_argument('--end', type=float, default=600.)
    run(parser.parse_args())
