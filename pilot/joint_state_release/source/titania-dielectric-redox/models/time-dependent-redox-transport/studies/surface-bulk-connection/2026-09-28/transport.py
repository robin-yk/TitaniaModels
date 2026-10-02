"""Apply the existing conservative transport solver to the joint-state pilot."""
import argparse
import hashlib
import importlib.util
import json
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import connection as c

OLD = c.PROJECT / 'models/vacancy-transport/time-dependent-redox-transport/versions/0.1.0'
SOURCE = c.PROJECT / 'models/vacancy-distribution/atomic-layer-and-bulk-defect-distribution/studies/site-population/reduction-source.csv'


def source():
    data = np.genfromtxt(SOURCE, delimiter=',', names=True)
    t, r = data['time_s'], data['rate_umol_g_s']
    te, re = np.r_[0., t, 1800.], np.r_[r[0], r, r[-1]]
    cum = np.r_[0., np.cumsum(np.diff(te) * (re[1:] + re[:-1]) / 2)]

    def amount(time):
        j = min(np.searchsorted(te, time, side='right') - 1, len(te) - 2)
        dt = time - te[j]
        return cum[j] + re[j] * dt + .5 * (re[j+1] - re[j]) / (te[j+1] - te[j]) * dt**2

    return amount, dict(recorded_start_s=float(t[0]), recorded_end_s=float(t[-1]),
                        recorded_integral_umol_g=float(np.trapezoid(r, t)),
                        extended_integral_umol_g=float(amount(1800.)),
                        same_cycle_reported_umol_g=94.8, nominal_distribution_umol_g=94.,
                        first_rate=float(r[0]), last_rate=float(r[-1]),
                        time_s=te.tolist(), rate_umol_g_s=re.tolist())


def run(args):
    adapter = ModuleType('pilot')
    adapter.HERE, adapter.T, adapter.pt = OLD, c.T, c.pt
    adapter.build = lambda n: c.build(n, args.delta, args.binding)
    adapter.quantities, adapter.source = c.quantities, source
    sys.modules['pilot'] = adapter
    spec = importlib.util.spec_from_file_location('original_transport', OLD / 'transient.py')
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    D = 1e13 * .3**2 / 6 * np.exp(-args.barrier / (c.pt.KB * c.T))
    output = c.HERE / f'transport-{args.barrier:.2f}-n{args.bulk}-dt{args.max_step:g}.json'
    engine.run(SimpleNamespace(bulk=args.bulk, D=D, end=1800., max_step=args.max_step,
                               relax=True, null_tail=True, output=str(output)))
    result = json.loads(output.read_text())
    result['model'] = 'Joint surface vacancy/electron states with ideal free-polaron bulk'
    result['assumed_parameters'] = dict(delta_eV=args.delta, binding_eV=args.binding,
                                       migration_barrier_eV=args.barrier,
                                       attempt_frequency_s=1e13, hop_length_nm=.3)
    for path in (c.HERE / 'connection.py', c.HERE / 'transport.py', SOURCE):
        result['code_sha256'][str(path.relative_to(c.PROJECT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    m, lay = c.build(args.bulk, args.delta, args.binding)
    sol = m.solve(result['source']['extended_integral_umol_g'], c.T)
    result['same_inventory_equilibrium'] = c.quantities(m, lay, sol.x, sol.phi)
    checks = result['checks']
    assert checks['max_vacancy_balance_umol_g'] < 1e-7
    assert checks['max_electron_balance_umol_g'] < 1e-7
    assert checks['min_transport_dissipation'] > -1e-8
    assert max(checks['closed_free_energy_increments']) < 1e-6
    output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--barrier', type=float, required=True)
    parser.add_argument('--bulk', type=int, default=80)
    parser.add_argument('--max-step', type=float, default=30.)
    parser.add_argument('--delta', type=float, default=-1.2)
    parser.add_argument('--binding', type=float, default=.8)
    run(parser.parse_args())
