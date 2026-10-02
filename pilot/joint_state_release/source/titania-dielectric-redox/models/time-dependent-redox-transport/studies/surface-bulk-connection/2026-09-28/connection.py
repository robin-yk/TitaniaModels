"""Joint surface vacancy/electron states connected to a finite R600 particle.

Isolated sensitivity pilot. Delta and binding are assumed bare local energies.
The archived free-energy engine supplies state tables and spherical Poisson
physics; a general covariance Hessian permits a nonadjacent electron plane.
"""

from pathlib import Path as _TitaniaPath
_TITANIA_PROJECT_ROOT = next(p for p in _TitaniaPath(__file__).resolve().parents if p.name == 'titania-dielectric-redox')
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import sys

import numpy as np
from scipy.linalg import solve
from scipy.special import logsumexp

HERE = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28'
PROJECT = _TITANIA_PROJECT_ROOT
ARCHIVE = _TITANIA_PROJECT_ROOT / 'models/source-snapshots/repository-baseline-9049a0a/pilot/titania-super-multiscale'
sys.path.insert(0, str(ARCHIVE))
from tofrange.engine import Model, COULOMB
from tofrange import particle as pt

T = 873.15
NOMINAL = 94.0
DIAMETER = 900.0
PLANES = {'BRI': 0, 'IPL': 1, 'SBR': 2}
MC_CACHE = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/versions/0.2.0/inputs/lnq_10dfef869565532e.json'


class ConnectedModel(Model):
    def _setup_electrostatics(self, domains, es):
        order = sorted(es['radii'], key=lambda s: -es['radii'][s])
        r = np.array([es['radii'][s] for s in order])
        assert np.all(np.diff(r) < 0)
        pos = {s: k for k, s in enumerate(order)}
        self.shell = np.array([pos[d['shell']] for d in domains])
        self.shell2 = np.array([pos[d.get('shell2', d['shell'])] for d in domains])
        c = 1 / np.diff(np.r_[0., 1 / r])
        self.es = dict(r=r, order=order, s=1 / r, L=len(r),
                       K=COULOMB / es['eps'] * 1e-6 * pt.NA * es['mass_g'],
                       main=c + np.r_[c[1:], 0.], off=-c[1:],
                       region_of_shell=np.zeros(len(r), dtype=int))

    def evaluate(self, Y, N, T, hessian=True):
        L = self.es['L']
        kT = pt.KB * T
        y, u = Y[:2], Y[2:]
        A = (self.logg - self.E / kT + y[0] * self.v + y[1] * self.e
             - u[self.shell, None] * self.zA - u[self.shell2, None] * self.zB)
        Z = logsumexp(A, axis=1)
        p = np.exp(A - Z[:, None])
        x = self.C[:, None] * p
        q = (np.bincount(self.shell, (x * self.zA).sum(1), minlength=L)
             + np.bincount(self.shell2, (x * self.zB).sum(1), minlength=L))
        Pu = self._Pmul(u, kT)
        g = np.r_[(x * self.v).sum() - N, (x * self.e).sum() - 2 * N, Pu - q]
        psi = self.C @ Z - N * y[0] - 2 * N * y[1] + .5 * u @ Pu
        if not hessian:
            return psi, x, g, q
        features = (self.v, self.e, -self.zA, -self.zB)
        ids = (np.zeros(len(self.C), int), np.ones(len(self.C), int),
               2 + self.shell, 2 + self.shell2)
        means = [(p * f).sum(1) for f in features]
        H = np.zeros((L + 2, L + 2))
        for a in range(4):
            for b in range(4):
                cov = self.C * ((p * features[a] * features[b]).sum(1) - means[a] * means[b])
                np.add.at(H, (ids[a], ids[b]), cov)
        j = np.arange(L) + 2
        H[j, j] += kT / self.es['K'] * self.es['main']
        H[j[:-1], j[1:]] += kT / self.es['K'] * self.es['off']
        H[j[1:], j[:-1]] += kT / self.es['K'] * self.es['off']
        return psi, x, g, q, H

    def _solve_global(self, N, T, tol, fixed):
        if fixed:
            raise ValueError('Frozen populations are outside this pilot')
        _, yL, chi, _, _, _ = self._solve_local(N, T, tol, {})
        Y = np.r_[yL + 2 * chi[0], -chi[0], np.zeros(self.es['L'])]
        threshold = max(tol * N, 2e-10)
        for it in range(1, 201):
            psi, x, g, q, H = self.evaluate(Y, N, T)
            if np.max(np.abs(g)) < threshold:
                break
            dY = solve(H, -g, assume_a='pos')
            slope = g @ dY
            step = 1.
            while step > 1e-12:
                newpsi, _, newg, _ = self.evaluate(Y + step * dY, N, T, False)
                if (newpsi <= psi + 1e-4 * step * slope + 1e-12
                        or np.max(np.abs(newg)) < .8 * np.max(np.abs(g))):
                    break
                step *= .5
            Y += step * dY
        else:
            raise RuntimeError(f'Equilibrium failed: {np.max(np.abs(g))}')
        self.last_Y = Y.copy()
        phi = Y[2:] * pt.KB * T
        poisson = float(np.max(np.abs(self.potentials(q) - phi)))
        return x, Y[:2], Y[2:], it, float(g[0]), float(g[1]), q, poisson


def joint_states(delta, binding, electron_energy=0., method='PBE', single=False):
    """Abstract cell: one O capacity plus two effective Ti electron capacities.

    A neutral defect costs delta relative to an isolated bulk vacancy and two
    bulk Ti electrons. Removing one electron loses binding/2. g=2 counts
    effective position choices, not spin. Three neutral representatives are
    an explicit coarse state-count assumption, not a lattice tiling.
    """
    gaps = [0.] if single else ([0., .28, .23] if method == 'PBE' else [0., .28, .27])
    v = np.array([0, 0, 0, 1, 1] + [1] * len(gaps))
    e = np.array([0, 1, 2, 0, 1] + [2] * len(gaps))
    E = np.array([0., electron_energy, 2 * electron_energy,
                  delta + binding - 2 * electron_energy,
                  delta + .5 * binding - electron_energy] + [delta + x for x in gaps])
    return dict(E=E.tolist(), v=v.tolist(), e=e.tolist(), g=[1, 2, 1, 1, 2] + [1] * len(gaps))


def build(n_bulk=80, delta=-1.2, binding=.8, bulk='ideal', eps=64., method='PBE',
          single=False, gauge=0.):
    R = DIAMETER / 2
    f110 = .75
    doms, radii, caps = [], {}, []
    rawO, rawTi = pt.layer_sites(DIAMETER, 4)
    oxy = [(k, s, c * f110, z) for k, s, c, z in rawO]
    ti = np.array(rawTi) * f110
    Cjoint = min(oxy[0][2], ti[1] / 2)

    def add(d, depth, O=0., Ti=0.):
        if d['C'] < 1e-12:
            return
        d['region'] = 0
        v = np.asarray(d['v'])
        e = np.asarray(d.get('e', np.zeros(len(v))))
        d['E'] = (np.asarray(d['E']) + gauge * (e - 2 * v)).tolist()
        radii[d['shell']] = R - depth
        doms.append(d)
        caps.append((O, Ti))

    def two(C, E, sh, tag, z, electron=False):
        add(dict(C=C, E=[0., E], v=[0, 0] if electron else [0, 1],
                 e=[0, 1] if electron else [0, 0], shell=sh, tag=tag), z,
            O=0. if electron else C, Ti=C if electron else 0.)

    for k, s, c, z in oxy:
        sh = (k, PLANES[s])
        if k == 1 and s == 'BRI':
            joint_id = len(doms)
            add(dict(C=Cjoint, shell=sh, shell2=(2, 1), tag='joint',
                     **joint_states(delta, binding, method=method, single=single)), z,
                O=Cjoint, Ti=2 * Cjoint)
            two(c - Cjoint, delta + binding, sh, 'bridging', z)
        else:
            two(c, 0., sh, (k, s), z)
    # Restore all L1 Ti. Reserve two electron capacities per joint cell at L2.
    for k, c in enumerate(ti, 1):
        two(c - (2 * Cjoint if k == 2 else 0.), 0., (k, 1), 'Ti',
            (k - 1) * pt.D110 + pt.H_BRI, electron=True)
    top = 4 * pt.D110
    edges = top + np.r_[0., np.geomspace(.02, R - top, n_bulk)]
    edges[-1] = R
    orem = pt.O_TOTAL - sum(c for _, _, c, _ in oxy)
    tirem = pt.TI_TOTAL - sum(ti)
    if bulk == 'historical':
        cache = json.loads(MC_CACHE.read_text())
        lnq = np.asarray(cache['lnQ'])
    for j, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        frac = pt.shell_fraction(R, lo, hi) / pt.shell_fraction(R, top, R)
        z, sh = .5 * (lo + hi), ('bulk', j)
        if bulk == 'ideal':
            two(orem * frac, 0., sh, 'bulk', z)
        else:
            add(dict(C=orem * frac / 700, E=(-pt.KB * T * lnq).tolist(),
                     v=np.arange(len(lnq)).tolist(), shell=sh, tag='bulk'), z, O=orem * frac)
        two(tirem * frac, 0., sh, 'Ti', z, electron=True)
    m = ConnectedModel(doms, dict(radii=radii, eps=eps, mass_g=pt.particle_mass(DIAMETER)))
    Ocap = np.zeros(m.es['L'])
    Tcap = np.zeros_like(Ocap)
    for d, (o, e) in enumerate(caps):
        Ocap[m.shell[d]] += o
        Tcap[m.shell2[d]] += e
    lay = SimpleNamespace(Ocap=Ocap, Tcap=Tcap, joint_id=joint_id, edges=edges,
                          delta=delta, binding=binding, bulk=bulk, eps=eps, method=method,
                          Cjoint=Cjoint, c_bri=oxy[0][2])
    assert np.isclose(Ocap.sum(), pt.O_TOTAL, rtol=1e-14)
    assert np.isclose(Tcap.sum(), pt.TI_TOTAL, rtol=1e-14)
    return m, lay


def quantities(m, lay, x, phi):
    L = m.es['L']
    vac = np.bincount(m.shell, (x * m.v).sum(1), minlength=L)
    el = np.bincount(m.shell2, (x * m.e).sum(1), minlength=L)
    pools = {'bridging': 0., 'other_explicit': 0., 'bulk': 0.}
    for d, tag in enumerate(m.tags):
        if tag == 'Ti':
            continue
        key = 'bulk' if tag == 'bulk' else 'bridging' if tag in ('joint', 'bridging') else 'other_explicit'
        pools[key] += float(x[d] @ m.v[d])
    jd = lay.joint_id
    jv = float(x[jd] @ m.v[jd])
    neutral = float(x[jd][(m.v[jd] == 1) & (m.e[jd] == 2)].sum())
    return dict(depth_nm=(DIAMETER / 2 - m.es['r']).tolist(), vacancy_umol_g=vac.tolist(),
                electron_umol_g=el.tolist(), oxygen_capacity_umol_g=lay.Ocap.tolist(),
                titanium_capacity_umol_g=lay.Tcap.tolist(),
                vacancy_pct=(100 * vac / lay.Ocap).tolist(),
                Ti3_pct=np.divide(100 * el, lay.Tcap, out=np.zeros(L), where=lay.Tcap > 0).tolist(),
                phi_eV_per_e=np.asarray(phi).tolist(), theta_BRI_pct=100 * pools['bridging'] / lay.c_bri,
                reconstruction_fraction=0., pools_umol_g=pools,
                vacancy_total_umol_g=float(vac.sum()), electron_total_umol_g=float(el.sum()),
                joint_neutral_fraction=neutral / jv if jv > 0 else 0.,
                joint_state_probability=(x[jd] / m.C[jd]).tolist())


def equilibrium(delta, binding, n_bulk=80, bulk='ideal', **kwargs):
    m, lay = build(n_bulk, delta, binding, bulk, **kwargs)
    sol = m.solve(NOMINAL, T)
    q = quantities(m, lay, sol.x, sol.phi)
    return dict(delta_eV=delta, binding_eV=binding, bulk=bulk, n_bulk=n_bulk,
                mu_v_eV=sol.mu_eV, mu_e_eV=sol.mu_e, iterations=sol.iterations,
                mass_error=sol.inventory_residual, electron_error=sol.charge_residual,
                poisson_error_eV=sol.poisson_residual, **q), m, lay, sol


def validate():
    out = {}
    r, m, lay, sol = equilibrium(-1.2, .8)
    out['capacity_O_error'] = float(lay.Ocap.sum() - pt.O_TOTAL)
    out['capacity_Ti_error'] = float(lay.Tcap.sum() - pt.TI_TOTAL)
    out['mass_error'] = r['mass_error']
    out['electron_error'] = r['electron_error']
    out['poisson_error_eV'] = r['poisson_error_eV']
    _, _, grad, _, H = m.evaluate(m.last_Y, NOMINAL, T)
    rng = np.random.default_rng(2819)
    d = rng.normal(size=len(grad)); d /= np.linalg.norm(d)
    h = 1e-5
    gp = m.evaluate(m.last_Y + h * d, NOMINAL, T, False)[2]
    gm = m.evaluate(m.last_Y - h * d, NOMINAL, T, False)[2]
    out['hessian_directional_rel_error'] = float(np.linalg.norm((gp - gm) / (2 * h) - H @ d) / np.linalg.norm(H @ d))
    out['hessian_symmetry_error'] = float(np.max(np.abs(H - H.T)))
    rg, mg, lg, sg = equilibrium(-1.2, .8, gauge=.37)
    out['gauge_population_max_error'] = float(np.max(np.abs(sg.x - sol.x)))
    out['gauge_mu_v_error'] = float(sg.mu_eV - sol.mu_eV + .74)
    out['gauge_mu_e_error'] = float(sg.mu_e - sol.mu_e - .37)
    probs = np.array(r['joint_state_probability'][5:8])
    expected = np.exp(-np.array([0., .28, .23]) / (pt.KB * T)); expected /= expected.sum()
    out['DFT_conditional_probability_error'] = float(np.max(np.abs(probs / probs.sum() - expected)))
    # Independent-site limit of the cell requires ONE neutral representative.
    st = joint_states(0., 0., single=True)
    z = np.exp(-np.asarray(st['E']) / (pt.KB * T) + .2 * np.array(st['v']) - .7 * np.array(st['e'])) * st['g']
    out['independent_cell_partition_error'] = float(abs(z.sum() - (1 + np.exp(.2)) * (1 + np.exp(-.7)) ** 2))
    qtest = np.zeros(m.es['L']); qtest[0] = 1.; qtest[4] = -1.
    phi = m.potentials(qtest)
    exact = m.es['K'] * (1 / m.es['r'][4] - 1 / m.es['r'][0])
    out['capacitor_energy_error'] = float(abs(.5 * qtest @ phi - .5 * exact))
    r240, _, _, _ = equilibrium(-1.2, .8, n_bulk=240)
    out['grid_80_240_BRI_pct_difference'] = abs(r240['theta_BRI_pct'] - r['theta_BRI_pct'])
    out['grid_80_240_bulk_umol_difference'] = abs(r240['pools_umol_g']['bulk'] - r['pools_umol_g']['bulk'])
    lnq = np.array(json.loads(MC_CACHE.read_text())['lnQ'])
    out['MC_single_vacancy_reference_error_eV'] = float(-pt.KB * T * lnq[1] + pt.KB * T * np.log(700))
    checks = {k: abs(v) < (1e-5 if k.startswith('grid') else 1e-7) for k, v in out.items()}
    # Grid differences are measured in percentage points / umol g-1.
    checks['grid_80_240_BRI_pct_difference'] = out['grid_80_240_BRI_pct_difference'] < .05
    checks['grid_80_240_bulk_umol_difference'] = out['grid_80_240_bulk_umol_difference'] < .05
    return dict(values=out, checks=checks, passed=all(checks.values()))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scan', action='store_true')
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()
    if args.test:
        tests = validate()
        (_TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/validation.json').write_text(json.dumps(tests, indent=2) + '\n')
        print(json.dumps(tests, indent=2), flush=True)
        if not tests['passed']:
            raise SystemExit('Validation failure')
    if args.scan:
        rows = []
        for delta in (-.4, -.8, -1.2, -1.6):
            for binding in (0., .4, .8):
                r, _, _, _ = equilibrium(delta, binding)
                rows.append(r)
                print(delta, binding, 'BRI%', round(r['theta_BRI_pct'], 3),
                      'bulk umol/g', round(r['pools_umol_g']['bulk'], 3),
                      'neutral%', round(100 * r['joint_neutral_fraction'], 2), flush=True)
        for bulk in ('historical',):
            r, _, _, _ = equilibrium(-1.2, .8, bulk=bulk)
            rows.append(r)
        for kwargs in (dict(eps=107.), dict(method='HSE')):
            r, _, _, _ = equilibrium(-1.2, .8, **kwargs)
            r['sensitivity'] = kwargs
            rows.append(r)
        payload = dict(model='Surface-bulk joint-state connection pilot 0.1.0',
                       sample='R600', T_K=T, diameter_nm=DIAMETER, inventory_umol_g=NOMINAL,
                       assumed_bare_energies=True, rows=rows,
                       sha256={str(p.relative_to(PROJECT)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__), _TITANIA_PROJECT_ROOT / 'models/source-snapshots/repository-baseline-9049a0a/pilot/titania-super-multiscale/tofrange/engine.py',
                                         _TITANIA_PROJECT_ROOT / 'models/source-snapshots/repository-baseline-9049a0a/pilot/titania-super-multiscale/tofrange/particle.py', MC_CACHE)})
        (_TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/equilibrium.json').write_text(json.dumps(payload, indent=2) + '\n')


if __name__ == '__main__':
    main()
