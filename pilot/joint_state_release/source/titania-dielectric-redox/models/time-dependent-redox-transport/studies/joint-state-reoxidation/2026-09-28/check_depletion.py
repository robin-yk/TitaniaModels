"""Test the surface nonnegativity boundary at a saved reoxidation stop.

Interior vacancy inventories are held fixed. Surface activity is reduced, while
electrons and potential re-equilibrate under total charge neutrality. This tests
the instantaneous limiting supply of that interior profile, not a global maximum.
"""
import argparse
import json
import warnings

import numpy as np
from scipy.linalg import solve, lstsq, LinAlgWarning

import reoxidation as r


def constrained(tr, start, interior, surface_alpha):
    L = tr.L
    Y = start.copy()
    Y[0] = surface_alpha

    def system(Y, jacobian):
        out = tr.state(Y, jacobian)
        g = out[0]
        f = np.r_[Y[0]-surface_alpha, g[1:L]-interior,
                  g[L]-2*g[:L].sum(), g[L+1:]]
        if not jacobian:
            return f
        H = out[3]
        unit = np.zeros(tr.size)
        unit[0] = 1.
        return f, np.vstack([unit, H[1:L], H[L]-2*H[:L].sum(0), H[L+1:]])

    for it in range(50):
        f, A = system(Y, True)
        if np.max(np.abs(f)) < 1e-10:
            return Y, it, float(np.max(np.abs(f)))
        rows = np.maximum(np.max(np.abs(A), axis=1), 1e-16)
        with warnings.catch_warnings():
            warnings.simplefilter('error', LinAlgWarning)
            try:
                d = solve(A/rows[:, None], -f/rows, check_finite=False)
            except (LinAlgWarning, np.linalg.LinAlgError):
                d = lstsq(A/rows[:, None], -f/rows, cond=1e-12,
                          check_finite=False, lapack_driver='gelsd')[0]
        length = min(1., 20/max(np.max(np.abs(d)), 1e-30))
        while length > 1e-8:
            Yn = Y + length*d
            if np.max(np.abs(system(Yn, False))) < np.max(np.abs(f)):
                Y = Yn
                break
            length *= .5
        else:
            raise RuntimeError('Constrained surface-depletion solve failed')
    raise RuntimeError('Constrained surface-depletion iteration limit')


def main(args):
    saved = json.loads((r.HERE/args.input).read_text())
    p = saved['parameters']
    e = r.engine(p['delta_eV'], p['binding_eV'])
    tr = e.Transport(saved['n_bulk'], saved['D_nm2_s'])
    last = saved['checkpoints'][-1]
    original = np.array(last['full_state_Y'])
    n = tr.state(original, False)[0][:tr.L]
    demand = r.measured()[1](saved['end_time_s'])
    results = []
    Y = original.copy()
    for reduction in [0., 10., 20., 30.]:
        Y, iterations, residual = constrained(tr, Y, n[1:], original[0]-reduction)
        g = tr.state(Y, False)[0]
        theta = g[:tr.L]/tr.cap
        eta = Y[:tr.L]-np.log(theta/(1-theta))
        inward = tr.flux(Y, g[:tr.L])[0]
        # J_out = G B(eta_0-eta_1) theta_1 when theta_0 -> 0.
        supply_limit = float(tr.G[0]*e.bernoulli(eta[0]-eta[1])[0]*theta[1])
        results.append(dict(surface_alpha_reduction=reduction,
            surface_vacancies_umol_g=float(g[0]), theta_surface=float(theta[0]),
            outward_flux_umol_g_s=float(-inward),
            vanishing_surface_limit_umol_g_s=supply_limit,
            imposed_demand_umol_g_s=demand,
            surface_inventory_derivative_umol_g_s=float(-inward-demand),
            max_constraint_residual=residual, iterations=iterations))
    assert results[-1]['surface_vacancies_umol_g'] < 1e-18
    assert abs(results[-1]['outward_flux_umol_g_s']-results[-2]['outward_flux_umol_g_s']) < 1e-8
    assert abs(results[-1]['outward_flux_umol_g_s']-results[-1]['vanishing_surface_limit_umol_g_s']) < 1e-8
    assert results[-1]['surface_inventory_derivative_umol_g_s'] < 0
    out = dict(input=args.input, stop_time_s=saved['end_time_s'],
        interpretation='For this saved interior profile and unchanged mobility, the imposed demand '
        'points outside the nonnegative surface-inventory domain. This diagnoses a supply limit '
        'at the observed stop; it does not calculate CO evolution after depletion.', cases=results)
    (r.HERE/args.output).write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='reox-2.00-n240-dt0.05.json')
    parser.add_argument('--output', default='depletion-check.json')
    main(parser.parse_args())
