"""The vacancy-distribution workspace reads paper_outputs/tof_range.json
and draws it. The gates: that file is what the one-model run wrote, it is
current, the page carries exactly it, and the surface panel places
vacancies from its values by the same rule in Python and in the browser."""

import csv
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
TOF = ROOT / 'paper_outputs' / 'tof_range.json'
PILOT = ROOT / 'pilot' / 'titania-super-multiscale' / 'outputs'
PAGE = ROOT / 'docs' / 'index.html'


def test_latest_release_data_is_embedded_without_changes():
    import re
    html = PAGE.read_text()
    match = re.search(r'<script id="joint-state-data" type="application/json">(.*?)</script>', html, re.S)
    assert match
    stored = json.loads((ROOT / 'pilot/joint_state_release/page_data.json').read_text())
    assert json.loads(match.group(1)) == stored
    assert stored['diameter_nm'] == 900
    assert stored['temperature_K'] == 873.15
    for panel in ['jrA', 'jrB', 'jrC', 'jrD']:
        assert f'id="{panel}"' in html
    assert 'window.TimeDependentDistribution' in html
    assert 'window.VacancyDistribution=' not in html


def test_latest_release_reoxidation_fills_vacancies():
    data = json.loads((ROOT / 'pilot/joint_state_release/page_data.json').read_text())
    ox = data['reoxidation']
    assert ox[0]['vacancy_total_umol_g'] == data['reduction'][-1]['vacancy_total_umol_g']
    assert all(b['vacancy_total_umol_g'] <= a['vacancy_total_umol_g'] + 1e-9 for a, b in zip(ox, ox[1:]))
    assert abs(ox[-1]['vacancy_total_umol_g'] - 14.19409) < 1e-4
    html = PAGE.read_text()
    assert 'oxygen incorporated into the solid' in html
    assert 'cumulative oxygen removal as CO' not in html


def test_release_site_fractions_and_charge_balance():
    import math
    data = json.loads((ROOT / 'pilot/joint_state_release/page_data.json').read_text())
    for history in ('reduction', 'reoxidation'):
        for row in data[history]:
            assert row['electron_total_umol_g'] == pytest.approx(2 * row['vacancy_total_umol_g'], abs=1e-6)
            if 'depth_nm' not in row:
                continue
            for fraction, amount, capacity in [('Ti3_pct', 'electron_umol_g', 'titanium_capacity_umol_g'), ('vacancy_pct', 'vacancy_umol_g', 'oxygen_capacity_umol_g')]:
                assert len(row[fraction]) == len(row[amount]) == len(row[capacity]) == len(row['depth_nm'])
                for y, n, c in zip(row[fraction], row[amount], row[capacity]):
                    assert all(math.isfinite(v) for v in (y, n, c))
                    assert -1e-10 <= n <= c + 1e-8
                    if c > 0:
                        assert y == pytest.approx(100 * n / c, abs=1e-10)
                    else:
                        assert n == pytest.approx(0, abs=1e-12)


def test_release_plot_does_not_invent_log_floor():
    source = (ROOT / 'web/distribution_release.js').read_text()
    assert 'Math.max(1e-10,y)' not in source
    assert 's.capacity[i]>0' in source
    assert "key==='Ti3_pct'?'titanium_capacity_umol_g':'oxygen_capacity_umol_g'" in source


@pytest.fixture(scope='module')
def doc():
    return json.loads(TOF.read_text())


def test_regenerating_the_json_is_byte_identical():
    before = TOF.read_bytes()
    r = subprocess.run([sys.executable, 'scripts/export_tof_range.py'],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert TOF.read_bytes() == before, 'run scripts/export_tof_range.py and commit'


def test_the_ranges_are_the_run_summary_as_written(doc):
    with (PILOT / 'sample_tof_range.csv').open() as fh:
        rows = {r['sample']: r for r in csv.DictReader(fh)}
    for s in doc['samples']:
        assert s['summary'] == rows[s['sample']]


def test_every_case_is_the_run_row_to_four_figures(doc):
    grid = doc['axes']
    eps = {f"{e['value']:g} ({e['key']})": i for i, e in enumerate(grid['eps'])}
    by = {s['sample']: s for s in doc['samples']}
    with (PILOT / 'model_cases.csv').open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 6 * 720 * 5
    for r in rows[::97]:
        i = grid['energy_map'].index(r['energy_map'])
        i = i * len(grid['cutoff_nm']) + grid['cutoff_nm'].index(float(r['cutoff_nm']))
        i = i * len(grid['dG_eV']) + grid['dG_eV'].index(float(r['dG_eV']))
        i = i * len(grid['f110']) + grid['f110'].index(float(r['f110']))
        i = i * len(grid['eps']) + eps[r['eps']]
        c, j = by[r['sample']]['cases'], doc['reactive'].index(r['reactive'])
        assert c['sites'][i][j] == pytest.approx(float(r['N_react_umol_g']), rel=5e-4)
        assert c['theta'][i] == pytest.approx(float(r['theta_bri']), rel=5e-4)
        assert c['below'][i][j] == (r['site_class'] != 'OK')
        tof = by[r['sample']]['rate_co_umol_g_s'] / c['sites'][i][j]
        assert tof == pytest.approx(float(r['TOF_s_1']), rel=1e-3)


def test_the_page_carries_the_json_unchanged(doc):
    assert json.dumps(doc, separators=(',', ':')) in PAGE.read_text(), \
        'the page is stale - run scripts/build_site.py'


RENDER = ROOT / 'web' / 'render'
HARNESS = ROOT / 'tests' / 'js' / 'slab_harness.js'
sys.path.insert(0, str(ROOT / 'scripts'))
import slab_atoms  # noqa: E402


def _fractions(doc, s, i, f110):
    cap = doc['capacity_umol_g']['%g' % f110]
    pools = dict(zip(doc['pools'], s['cases']['pools'][i]))
    return dict(BRI=s['cases']['theta'][i], IPL=pools['basal_L1'] / cap['basal_L1'],
                SBR=pools['L1_subbridging'] / cap['L1_subbridging'],
                L24=pools['subsurface_L2_4'] / cap['subsurface_L2_4'])


def _f110(doc, i):
    ax = doc['axes']
    return ax['f110'][(i // len(ax['eps'])) % len(ax['f110'])]


def test_the_atom_list_regenerates_from_the_cif():
    before = (RENDER / 'slab_atoms.json').read_bytes()
    r = subprocess.run([sys.executable, 'scripts/slab_atoms.py'], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (RENDER / 'slab_atoms.json').read_bytes() == before, 'run scripts/slab_atoms.py and commit'


def test_every_o_and_ti_site_of_the_slab_is_counted_once():
    at = json.loads((RENDER / 'slab_atoms.json').read_text())
    nx, ny = at['cells']
    per = nx * ny * at['trilayers']
    from collections import Counter
    c = Counter(at['site'])
    assert c == {'BRI': nx * ny, 'IPL': 2 * nx * ny, 'SBR': nx * ny,
                 'L24': 4 * nx * ny * (at['trilayers'] - 1), 'Ti': 2 * per}
    assert len({tuple(x) for x in at['xyz']}) == len(at['xyz'])


def test_the_browser_places_the_same_vacancies(doc):
    """web/slab_canvas.js against scripts/slab_atoms.py: fractions from the
    stored case and the vacant sites, for every sample on a spread of
    parameter points."""
    import shutil
    assert shutil.which('node') is not None, 'node is required for this gate'
    npts = len(doc['samples'][0]['cases']['theta'])
    points = list(range(0, npts, 37)) + [npts - 1]
    r = subprocess.run(['node', str(HARNESS), str(ROOT), json.dumps(points)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    js = json.loads(r.stdout)
    at = json.loads((RENDER / 'slab_atoms.json').read_text())
    n = 0
    for s in doc['samples']:
        for i in points:
            f = _fractions(doc, s, i, _f110(doc, i))
            got = js['%s:%d' % (s['sample'], i)]
            assert got['fractions'] == pytest.approx(f, rel=1e-15), (s['sample'], i)
            assert got['vacant'] == slab_atoms.vacant(at, f), (s['sample'], i)
            n += len(got['vacant'])
    assert n > 0


def test_a_larger_fraction_keeps_every_vacancy_of_a_smaller_one(doc):
    """Fixed site ranks: when every class fraction of one case is at most
    that of another, its vacant sites are a subset of the other's."""
    at = json.loads((RENDER / 'slab_atoms.json').read_text())
    pairs = 0
    for i in (0, 123, 400):
        fr = [_fractions(doc, s, i, _f110(doc, i)) for s in doc['samples']]
        for f in fr:
            for g in fr:
                if f is not g and all(f[c] <= g[c] for c in f):
                    assert set(slab_atoms.vacant(at, f)) <= set(slab_atoms.vacant(at, g))
                    pairs += 1
    assert pairs > 0


def test_the_page_carries_the_renders():
    html = PAGE.read_text()
    assert 'RENDER:' not in html and '/*RENDERDATA*/' not in html and '/*SLAB*/' not in html
    for name in ('particle_default.webp', 'particle_mask.png', 'slab_atoms.json'):
        assert (RENDER / name).exists(), name
    assert html.count('data:image/webp;base64,') >= 4
    at = json.loads((RENDER / 'slab_atoms.json').read_text())
    assert json.dumps(at, separators=(',', ':')) in html
