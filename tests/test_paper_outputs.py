"""The published numbers and the package must not drift apart.

Everything the manuscript quotes is produced by one script, so the gates
are: the script is deterministic, its output is committed and current,
the page shows exactly that output, and the values in it are what a fresh
solve returns rather than a stale copy.
"""

import csv
import json
import pathlib
import subprocess
import sys

import pytest

from solidgas import activeset as A

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'paper_outputs'
PAGE = ROOT / 'docs' / 'index.html'
FEED = {'CO2': 10, 'H2': 10, 'N2': 80}


def _rows(name):
    with (OUT / name).open() as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope='module')
def site():
    with (OUT / 'site_data.json').open() as fh:
        return json.load(fh)


# ------------------------------------------------------------- it is current

def test_regenerating_reproduces_the_committed_outputs():
    """Byte for byte. A stale table is the failure mode worth a gate."""
    before = {p.name: p.read_bytes() for p in sorted(OUT.iterdir())}
    r = subprocess.run([sys.executable, 'scripts/reproduce_paper.py'],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    after = {p.name: p.read_bytes() for p in sorted(OUT.iterdir())}
    assert set(before) == set(after)
    stale = [k for k in before if before[k] != after[k]]
    assert not stale, ('run scripts/reproduce_paper.py and commit: %s' % stale)


def test_every_output_the_manuscript_needs_is_there():
    have = {p.name for p in OUT.iterdir()}
    assert {'equilibrium_conversion.csv', 'phase_stability.csv',
            'vacancy_population.csv', 'dielectric_correlation.csv',
            'site_data.json'} <= have


# ------------------------------------------------------- it is what it says

def test_the_conversion_table_is_a_fresh_solve():
    for row in _rows('equilibrium_conversion.csv')[::7]:
        T_C = float(row['T_C'])
        r = A.solve(FEED, T_C + 273.15)
        assert float(row['X_CO2_pct']) == pytest.approx(
            r['conversion_CO2_pct'], rel=1e-9), T_C
        assert row['stable_phases'] == '+'.join(r['active_condensed_phases'])


def test_the_manuscript_conversion_at_600C_is_the_number_in_the_paper():
    """The paper quotes an equilibrium limit of about 38 percent."""
    r = A.solve(FEED, 873.15)
    assert r['conversion_CO2_pct'] == pytest.approx(38.0, abs=0.1)


def test_the_phase_stability_table_is_a_fresh_closed_form():
    for row in _rows('phase_stability.csv')[::7]:
        b = A.reduction_boundary(float(row['T_C']) + 273.15)
        assert float(row['critical_y_CO2_mol_pct']) == pytest.approx(
            b['y_CO2'] * 100.0, rel=1e-9)


def test_rutile_survives_the_whole_reported_range():
    for row in _rows('equilibrium_conversion.csv'):
        assert row['stable_phases'] == 'TiO2', row['T_C']
        assert float(row['margin_kJ_per_mol_O']) > 0


def test_the_dielectric_row_is_carried_not_recomputed():
    """It is the manuscript's own fit; the repository must not re-derive it
    from a file it does not have."""
    row = _rows('dielectric_correlation.csv')[0]
    assert 'not recomputed' in row['source']
    assert int(row['n_points']) == 4


# ------------------------------------------------------------- the page

def test_the_page_shows_exactly_what_the_script_produced(site):
    html = PAGE.read_text()
    assert json.dumps(site, separators=(',', ':')) in html, \
        'the page is stale - run scripts/build_site.py'


def test_the_page_rebuilds_byte_for_byte():
    """The committed page is what build_site.py writes from the committed
    sources, template text included. The MathML depends on the latex2mathml
    version, so requirements.txt pins the one the page was built with."""
    sys.path.insert(0, str(ROOT / 'scripts'))
    import build_site
    assert build_site.render() == PAGE.read_text(), \
        ('docs/index.html does not rebuild byte for byte - run '
         'scripts/build_site.py with the pinned requirements and commit')


def test_the_page_makes_no_external_request():
    html = PAGE.read_text()
    assert '<script src=' not in html
    assert '<link ' not in html
    # links out to the repository and to reference DOIs are not requests
    assert 'href="http' not in html.replace('href="https://github.com', '').replace(
        'href="https://doi.org/', '')


def test_the_browser_engines_are_mirrors_under_a_parity_gate():
    """The page solves, and that is the point - it is the tool. What must
    not happen is a second engine drifting from the first unwatched, so
    every browser engine has to be paired with a gate that holds it to the
    Python module it mirrors."""
    engines = {'activeset.js': 'test_activeset_port.py',
               'population.js': 'test_population_port.py',
               'slab_canvas.js': 'test_distribution_page.py'}
    for js, gate in engines.items():
        assert (ROOT / 'web' / js).exists(), js
        assert (ROOT / 'tests' / gate).exists(), \
            '%s has no parity gate; %s must exist' % (js, gate)
    inlined = sorted(p.name for p in (ROOT / 'web').iterdir()
                     if p.suffix == '.js')
    assert inlined == ['activeset.js', 'distribution_release.js', 'distribution_ui.js', 'figkit.js',
                       'figures_distribution.js', 'figures_population.js',
                       'figures_thermo.js', 'population.js', 'site_ui.js',
                       'slab_canvas.js', 'thermo_ui.js', 'ti_solver_page.js'], inlined


def test_the_manuscript_values_are_read_not_recomputed(site):
    """The page may solve whatever the user types, but the numbers the
    paper quotes come out of the committed JSON, so the two cannot differ."""
    html = PAGE.read_text()
    assert json.dumps(site, separators=(',', ':')) in html


def test_the_page_states_its_inputs_and_limits():
    text = ' '.join(PAGE.read_text().split())
    assert 'n<sub>iso</sub> is inferred from the measured rates' in text
    assert 'effective parameters of the net isolated-vacancy loss rate' in text


def test_the_site_data_and_the_tables_agree(site):
    csv_rows = _rows('equilibrium_conversion.csv')
    assert len(site['equilibrium']['rows']) == len(csv_rows)
    for a, b in zip(site['equilibrium']['rows'], csv_rows):
        assert a['T_C'] == pytest.approx(float(b['T_C']))
        assert a['conversion_CO2_pct'] == pytest.approx(float(b['X_CO2_pct']))
