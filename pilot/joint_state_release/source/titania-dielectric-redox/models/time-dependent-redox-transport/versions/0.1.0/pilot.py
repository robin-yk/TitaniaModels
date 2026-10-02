"""R600 source-constrained pilot; production sources are imported read-only.

Run with the scipy environment and a TitaniaModels checkout argument.
The equilibrium calculation is a benchmark, never a transient solution.
"""

from pathlib import Path as _TitaniaPath
_TITANIA_PROJECT_ROOT = next(p for p in _TitaniaPath(__file__).resolve().parents if p.name == 'titania-dielectric-redox')

from pathlib import Path
import hashlib
import json
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

HERE = (_TITANIA_PROJECT_ROOT / 'models/vacancy-transport/time-dependent-redox-transport/versions/0.1.0')
REPO = Path(sys.argv[1] if len(sys.argv) > 1 else '/private/tmp/titania-dion-20260926')
sys.path.insert(0, str(REPO / 'pilot/titania-super-multiscale'))
from tofrange import model as md, particle as pt

T = 873.15
DIAMETER = 900.0
SOURCE = (_TITANIA_PROJECT_ROOT / 'models/vacancy-distribution/atomic-layer-and-bulk-defect-distribution/studies/site-population/reduction-source.csv')


def source():
    raw = np.genfromtxt(SOURCE, delimiter=',', names=True)
    t, r = raw['time_s'], raw['rate_umol_g_s']
    # Constant endpoint extension is declared separately from the measured interval.
    te, re = np.r_[0., t, 1800.], np.r_[r[0], r, r[-1]]
    cumulative = np.r_[0., np.cumsum(np.diff(te) * (re[1:] + re[:-1]) / 2)]

    def amount(time):
        j = min(np.searchsorted(te, time, side='right') - 1, len(te) - 2)
        dt = time - te[j]
        slope = (re[j + 1] - re[j]) / (te[j + 1] - te[j])
        return float(cumulative[j] + re[j] * dt + slope * dt * dt / 2)

    return amount, dict(recorded_start_s=float(t[0]), recorded_end_s=float(t[-1]),
                        recorded_integral_umol_g=float(np.trapezoid(r, t)),
                        extended_integral_umol_g=amount(1800.),
                        same_cycle_reported_umol_g=94.8,
                        first_rate_umol_g_s=float(r[0]), last_rate_umol_g_s=float(r[-1]),
                        time_s=te.tolist(), rate_umol_g_s=re.tolist())


def build(n_bulk=240):
    return md.particle('LI_SBR1', DIAMETER, cutoff=.28, dG=0., f110=.75,
                       eps_r=64., s0=.2, n_bulk=n_bulk)


def quantities(m, lay, x, phi):
    L = m.es['L']
    vacancy = np.bincount(m.shell, (x * m.v).sum(1), minlength=L)
    electron = np.bincount(m.shell2, (x * m.e).sum(1), minlength=L)
    ocap, tcap = np.zeros(L), np.zeros(L)
    for d, tag in enumerate(m.tags):
        if tag == 'bulk':
            ocap[m.shell[d]] += 700 * m.C[d]
        elif tag == 'Ti':
            tcap[m.shell[d]] += m.C[d]
        elif d == lay.idx['cell']:
            ocap[m.shell[d]] += 2 * m.C[d]
            tcap[m.shell2[d]] += 4 * m.C[d]
        else:
            ocap[m.shell[d]] += m.C[d]
    sol = SimpleNamespace(model=m, x=x)
    th, cap, rec = md.surface(sol, lay)
    z = DIAMETER / 2 - m.es['r']
    return dict(depth_nm=z.tolist(), vacancy_umol_g=vacancy.tolist(),
                electron_umol_g=electron.tolist(), oxygen_capacity_umol_g=ocap.tolist(),
                titanium_capacity_umol_g=tcap.tolist(),
                vacancy_pct=(100 * vacancy / ocap).tolist(),
                Ti3_pct=np.divide(100 * electron, tcap, out=np.zeros(L), where=tcap > 0).tolist(),
                phi_eV_per_e=phi.tolist(), theta_BRI_pct=float(100 * th),
                reconstruction_fraction=float(rec),
                pools_umol_g={k: float(v) for k, v in md.populations(sol, lay).items()},
                vacancy_total_umol_g=float(vacancy.sum()), electron_total_umol_g=float(electron.sum()))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    from scipy.optimize import brentq
    amount, src = source()
    m, lay = build()
    zero = quantities(m, lay, np.zeros_like(m.E), np.zeros(m.es['L']))
    caps = np.array(zero['oxygen_capacity_umol_g'])
    src['BRI_capacity_reached_s'] = brentq(lambda t: amount(t) - caps[0], 0., 1800.)
    out = dict(sample='CH600 first reduction, corresponding R600 condition',
               diameter_nm=DIAMETER, T_K=T, duration_s=1800., initial_vacancy_umol_g=0.,
               initial_electron_umol_g=0., source=src,
               assumptions=dict(energy_set='unchanged LI_SBR1', f110=.75, eps_r=64.,
                                reconstruction_dG_eV=0., MC_cutoff_nm=.28,
                                H2_sccm=2.5, total_flow_sccm=50., sample_g=1.15),
               capacity=dict(BRI_umol_g=float(caps[0]), first_trilayer_O_umol_g=float(caps[:3].sum()),
                             first_four_trilayer_O_umol_g=float(caps[:12].sum()),
                             first_four_trilayer_Ti_umol_g=float(sum(zero['titanium_capacity_umol_g'][:12]))),
               equilibrium_benchmarks=[], repo_commit=subprocess.check_output(
                   ['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip(),
               source_sha256={str(SOURCE): sha(SOURCE)})
    for name in ('engine.py', 'model.py', 'particle.py', 'aggregates.py'):
        path = REPO / 'pilot/titania-super-multiscale/tofrange' / name
        out['source_sha256'][str(path)] = sha(path)
    for t in (60., 180., 600., 1800.):
        N = amount(t)
        sol = m.solve(N, T)
        row = quantities(m, lay, sol.x, sol.phi)
        row.update(time_s=t, target_umol_g=N, benchmark='instantaneous full equilibration')
        checks = dict(vacancy_error_umol_g=row['vacancy_total_umol_g'] - N,
                      electron_error_umol_g=row['electron_total_umol_g'] - 2 * N,
                      poisson_error_eV_per_e=sol.poisson_residual,
                      normalisation_error=float(np.max(np.abs(sol.x.sum(1) / m.C - 1))))
        assert abs(checks['vacancy_error_umol_g']) < 1e-6
        assert abs(checks['electron_error_umol_g']) < 1e-6
        assert checks['poisson_error_eV_per_e'] < 1e-5
        assert checks['normalisation_error'] < 1e-12
        row['checks'] = checks
        out['equilibrium_benchmarks'].append(row)
        print(json.dumps({k: row[k] for k in ('time_s', 'target_umol_g', 'theta_BRI_pct', 'pools_umol_g')}), flush=True)
    ((_TITANIA_PROJECT_ROOT / 'models/vacancy-transport/time-dependent-redox-transport/versions/0.1.0/benchmark.json')).write_text(json.dumps(out, indent=2) + '\n')
    for path, digest in out['source_sha256'].items():
        assert sha(path) == digest
    print('Source and equilibrium benchmark saved.', flush=True)


if __name__ == '__main__':
    main()
