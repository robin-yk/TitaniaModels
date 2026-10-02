"""Depth-dependent mobility pilot; archived free energy and flux law unchanged.

Same mobility in both treatment stages. Harmonic face averaging integrates the
piecewise-constant diffusion resistance across each cell-centre interval.
"""

from pathlib import Path as _TitaniaPath
_TITANIA_PROJECT_ROOT = next(p for p in _TitaniaPath(__file__).resolve().parents if p.name == 'titania-dielectric-redox')
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np

HERE = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/depth-barrier-test/2026-09-29'
PROJECT = _TITANIA_PROJECT_ROOT
CONNECTION = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28'
REOX = _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28'
sys.path[:0] = [str(REOX), str(CONNECTION)]
import connection as c
import transport as reduction
import reoxidation as reox
import check_depletion


def install(outer, inner, depth):
    original = reox.engine
    def factory(delta, binding):
        e = original(delta, binding)
        base = e.Transport
        class LayerTransport(base):
            def __init__(self, n_bulk=80, diffusivity_nm2_s=None):
                super().__init__(n_bulk, diffusivity_nm2_s)
                z = 450. - self.m.es['r']
                length = np.diff(z)
                outer_length = np.clip(depth-z[:-1], 0., length)
                Do, Di = 1e13*.3**2/6*np.exp(-np.array([outer, inner])/self.kT)
                face_D = length/(outer_length/Do+(length-outer_length)/Di)
                self.G *= face_D/self.D
                self.mobility = dict(outer_barrier_eV=outer, interior_barrier_eV=inner,
                    outer_depth_nm=depth, outer_D_nm2_s=float(Do), interior_D_nm2_s=float(Di),
                    face_D_nm2_s=face_D.tolist(), face_depth_nm=((z[:-1]+z[1:])/2).tolist())
        e.Transport = LayerTransport
        return e
    reox.engine = factory


def run(a):
    install(a.outer, a.inner, a.depth)
    tag = f'o{a.outer:g}-i{a.inner:g}-z{a.depth:g}-n{a.bulk}-r{a.reduction_step:g}-x{a.reox_step:g}'
    initial, final = [HERE/f'{s}-{tag}.json' for s in ('reduction', 'reox')]
    if initial.exists() or final.exists():
        raise FileExistsError(tag)
    e = reox.engine(-1.2, .8)
    e.source = reduction.source
    D = float(1e13*.3**2/6*np.exp(-a.outer/(c.pt.KB*c.T)))
    e.run(SimpleNamespace(bulk=a.bulk, D=D, end=1800., max_step=a.reduction_step,
        relax=True, null_tail=True, output=str(initial)))
    out = json.loads(initial.read_text())
    out['model'] = 'model-tio2-surface-bulk-joint-state v0.1.0 with depth-dependent mobility'
    out['assumed_parameters'] = dict(delta_eV=-1.2, binding_eV=.8,
        outer_barrier_eV=a.outer, interior_barrier_eV=a.inner, outer_depth_nm=a.depth,
        attempt_frequency_s=1e13, hop_length_nm=.3)
    out['mobility'] = e.Transport(a.bulk, D).mobility
    dependencies = [Path(__file__), _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/connection.py', _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/transport.py',
        _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28/reoxidation.py', _TITANIA_PROJECT_ROOT / 'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28/check_depletion.py', reox.OLD/'transient.py',
        reduction.SOURCE, reox.SOURCE]
    out['dependency_sha256'] = {str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in dependencies}
    checks = out['checks']
    assert checks['max_vacancy_balance_umol_g'] < 1e-7
    assert checks['max_electron_balance_umol_g'] < 1e-7
    assert checks['min_transport_dissipation'] > -1e-8
    assert max(checks['closed_free_energy_increments']) < 1e-6
    initial.write_text(json.dumps(out, indent=2)+'\n')
    reox.run(SimpleNamespace(initial=str(initial), output=str(final),
        max_step=a.reox_step, min_step=1e-7, end=600.))
    result = json.loads(final.read_text())
    result['mobility'] = out['mobility']
    result['dependency_sha256'] = out['dependency_sha256']
    final.write_text(json.dumps(result, indent=2)+'\n')
    if result['status'] != 'completed':
        check_depletion.main(SimpleNamespace(input=str(final), output=str(HERE/f'depletion-{tag}.json')))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--outer', type=float, default=1.6)
    p.add_argument('--inner', type=float, default=2.)
    p.add_argument('--depth', type=float, required=True)
    p.add_argument('--bulk', type=int, default=80)
    p.add_argument('--reduction-step', type=float, default=30.)
    p.add_argument('--reox-step', type=float, default=1.)
    run(p.parse_args())
