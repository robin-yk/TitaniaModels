"""Export the verified R600 study. Numerical arrays are copied without changes."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

DEST = Path(__file__).resolve().parent

def main(root):
    project = root / 'B2 working-paper/manuscripts-B2/titania-dielectric-redox'
    study = project / 'models/vacancy-transport/time-dependent-redox-transport/studies/depth-barrier/2026-09-29'
    reduction_path = study / 'reduction-o1-6-i2-2-z30-n32-119902.json'
    reox_path = study / 'reox-o1.6-i2.2-z30-n320-r2.5-x0.05.json'
    reduction = json.loads(reduction_path.read_text())
    reox = json.loads(reox_path.read_text())
    summary = json.loads((study / 'summary.json').read_text())
    source_root = DEST / 'source/titania-dielectric-redox'
    paths = [
        'models/time-dependent-redox-transport/studies/depth-barrier-test/2026-09-29/run_case.py',
        'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/connection.py',
        'models/time-dependent-redox-transport/studies/surface-bulk-connection/2026-09-28/transport.py',
        'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28/reoxidation.py',
        'models/time-dependent-redox-transport/studies/joint-state-reoxidation/2026-09-28/check_depletion.py',
        'models/time-dependent-redox-transport/versions/0.1.0/transient.py',
        'models/time-dependent-redox-transport/versions/0.1.0/pilot.py',
        'models/atomic-layer-and-bulk-defect-distribution/studies/site-population/reduction-source.csv',
        'models/atomic-layer-and-bulk-defect-distribution/studies/site-population/reoxidation-source.csv',
    ]
    archive = 'models/source-snapshots/repository-baseline-9049a0a/pilot/titania-super-multiscale/tofrange/'
    paths += [archive + name for name in ('__init__.py', 'engine.py', 'particle.py')]
    paths += ['models/time-dependent-redox-transport/versions/0.2.0/inputs/lnq_10dfef869565532e.json']
    provenance = {}
    def locate(relative):
        actual = relative.replace('models/time-dependent-redox-transport/', 'models/vacancy-transport/time-dependent-redox-transport/')
        actual = actual.replace('/studies/depth-barrier-test/', '/studies/depth-barrier/').replace('/studies/surface-bulk-connection/', '/studies/joint-state/')
        actual = actual.replace('models/atomic-layer-and-bulk-defect-distribution/', 'models/vacancy-distribution/atomic-layer-and-bulk-defect-distribution/')
        return project / actual
    for relative in paths:
        src = locate(relative)
        dst = source_root / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        provenance[relative] = hashlib.sha256(src.read_bytes()).hexdigest()
    # The reduction module uses the post-reorganization canonical paths.
    for old, new in [
        ('models/time-dependent-redox-transport/versions/0.1.0/transient.py', 'models/vacancy-transport/time-dependent-redox-transport/versions/0.1.0/transient.py'),
        ('models/time-dependent-redox-transport/versions/0.1.0/pilot.py', 'models/vacancy-transport/time-dependent-redox-transport/versions/0.1.0/pilot.py'),
        ('models/atomic-layer-and-bulk-defect-distribution/studies/site-population/reduction-source.csv', 'models/vacancy-distribution/atomic-layer-and-bulk-defect-distribution/studies/site-population/reduction-source.csv'),
    ]:
        dst = source_root / new
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(locate(old), dst)
    raw = DEST / 'results'
    raw.mkdir(exist_ok=True)
    for src in (reduction_path, reox_path, study/'summary.json'):
        shutil.copy2(src, raw/src.name)
    for name in ('depth-dependent-mobility.svg', 'depth-dependent-mobility.png'):
        shutil.copy2(study/name, DEST/name)
    fields = ['time_s', 'depth_nm', 'vacancy_pct', 'Ti3_pct', 'phi_eV_per_e',
              'vacancy_umol_g', 'electron_umol_g', 'oxygen_capacity_umol_g',
              'titanium_capacity_umol_g', 'theta_BRI_pct', 'vacancy_total_umol_g',
              'electron_total_umol_g', 'pools_umol_g']
    def points(run):
        rows = run['checkpoints']
        selected = set(range(len(rows))) if len(rows) <= 6 else {round(i*(len(rows)-1)/6) for i in range(7)}
        out = []
        for index, row in enumerate(rows):
            item = {k: row[k] for k in fields if k in row and not isinstance(row[k], list)}
            if index in selected:
                n = len(row['depth_nm'])
                keep = sorted(set(range(min(12, n))) | {round(i*(n-1)/107) for i in range(108)})
                for k in fields:
                    if k in row and isinstance(row[k], list) and len(row[k]) == n:
                        item[k] = [row[k][i] for i in keep]
            out.append(item)
        return out
    data = dict(schema_version=1, model=reduction['model'], temperature_K=reduction['T_K'],
                diameter_nm=reduction['diameter_nm'], parameters=reduction['assumed_parameters'],
                reduction=points(reduction), reoxidation=points(reox),
                reduction_source=reduction['source'], reoxidation_source=reox['source'],
                checks=dict(reduction=reduction['checks'], reoxidation=reox['checks']),
                cases=summary['cases'], source_sha256=provenance,
                scientific_scope='R600 conditional transport feasibility; measured rates are boundary inputs.')
    (DEST/'page_data.json').write_text(json.dumps(data, separators=(',', ':'))+'\n')
    (DEST/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('wiki_root', type=Path)
    main(parser.parse_args().wiki_root)
