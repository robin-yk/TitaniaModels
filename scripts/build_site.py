"""Assemble the site: one self-contained page, no external requests.

    python3 scripts/export_activeset.py   # the coefficient tables
    python3 scripts/oracle_tio.py         # the 80-digit reference
    python3 scripts/reproduce_paper.py    # the committed manuscript values
    python3 scripts/build_site.py         # -> docs/index.html

Both workspaces solve in the page. The equilibrium half runs the active-set
solver mirrored in web/activeset.js, and the population half integrates the
rate equation mirrored in web/population.js; parity gates hold each against
its Python original. The committed manuscript values ride along as
site_data.json so the figures and the paper cannot drift.
"""

import json
import re
import os
import base64
import xml.etree.ElementTree as ET
from site_math import typeset

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web')
DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(ROOT, 'docs', 'index.html')
SITE = os.path.join(ROOT, 'paper_outputs', 'site_data.json')
TOF = os.path.join(ROOT, 'paper_outputs', 'tof_range.json')

CASES_FEEDS = {
    'rwgs_1_1': {'CO2': 1, 'H2': 1},
    'h2_rich': {'CO2': 1, 'H2': 3},
    'product_mix': {'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2},
    'pure_h2': {'H2': 1},
    'pure_n2': {'N2': 1},
    'n2_10ppm_o2': {'N2': 999990, 'O2': 10},
}

PARTS = [
    ('/*DISTFIG*/', os.path.join(WEB, 'figures_distribution.js')),
    ('/*SLAB*/', os.path.join(WEB, 'slab_canvas.js')),
    ('/*DISTUI*/', os.path.join(WEB, 'distribution_release.js')),
    ('/*CSS*/', os.path.join(WEB, 'site.css')),
    ('/*FIGKIT*/', os.path.join(WEB, 'figkit.js')),
    ('/*ENGINE*/', os.path.join(WEB, 'activeset.js')),
    ('/*THERMO*/', os.path.join(WEB, 'figures_thermo.js')),
    ("/*POPENGINE*/", os.path.join(WEB, "population.js")),
    ("/*POPFIG*/", os.path.join(WEB, "figures_population.js")),
    ('/*PAGE*/', os.path.join(WEB, 'ti_solver_page.js')),
    ('/*THERMOUI*/', os.path.join(WEB, 'thermo_ui.js')),
    ('/*UI*/', os.path.join(WEB, 'site_ui.js')),
]


def read(path):
    with open(path) as fh:
        return fh.read()


def slim_reference():
    """Winner-level fields only: what the Validation tab diffs against."""
    doc = json.loads(read(os.path.join(DATA,
                                       'reference_results_high_precision.json')))
    rows = []
    for r in doc['rows']:
        if r['case'] not in CASES_FEEDS:
            continue
        rows.append({
            'case': r['case'], 'T_C': r['T_C'],
            'active': r['active_condensed_phases'],
            'gas_fractions': r['gas_fractions'],
            'r_per_mol_O': {s: rc['per_mol_O_kJ'] for s, rc
                            in r['inactive_phase_reduced_costs'].items()},
            'log10_traces': r['log10_trace_amounts'],
            'reduced_pct': r['reduced_pct'],
        })
    return {'precision': doc['precision'], 'method': doc['method'],
            'cases_feeds': CASES_FEEDS, 'rows': rows,
            'h2_h2o_boundary': doc['boundary_validation']['h2_h2o_boundary']}


RENDER = os.path.join(WEB, 'render')


def data_uri(name):
    kind = 'image/webp' if name.endswith('.webp') else 'image/png'
    with open(os.path.join(RENDER, name), 'rb') as fh:
        return 'data:%s;base64,%s' % (kind, base64.b64encode(fh.read()).decode())


def inline_renders(html):
    """Renders of the vacancy-distribution workspace. The grain image the
    page opens on goes straight into <img src>; the hover states, the hover
    mask and the (110) atom list go into one JSON block the script reads."""
    for name in re.findall(r'RENDER:([\w.]+)', html):
        html = html.replace('RENDER:' + name, data_uri(name))
    doc = dict(mask=data_uri('particle_mask.png'),
               particle={k: data_uri('particle_%s.webp' % k)
                         for k in ('default', 'surface', 'subsurface', 'bulk')},
               slab=json.loads(read(os.path.join(RENDER, 'slab_atoms.json'))))
    return html.replace('/*RENDERDATA*/', json.dumps(doc, separators=(',', ':')))


def render():
    """The page as a string; build() writes it and a gate compares it."""
    if not os.path.exists(SITE):
        raise SystemExit('run scripts/reproduce_paper.py first')
    html = typeset(read(os.path.join(WEB, 'template.html')))
    for name in ('gas-solid-equilibrium', 'vacancy-kinetics'):
        svg = read(os.path.join(WEB, 'schemes', name + '.svg'))
        root = ET.fromstring(svg)
        ns = {'s': 'http://www.w3.org/2000/svg'}
        # Workspace headings already identify the model. Remove the plate label,
        # repeated title and subtitle, retaining the physical drawing and labels.
        for parent in root.iter():
            for child in list(parent):
                if child.get('id') in ('text_1', 'text_2', 'text_3'):
                    parent.remove(child)
        root.set('viewBox', '0 55 360 225')
        root.set('height', '225pt')
        svg = ET.tostring(root, encoding='unicode')
        html = html.replace('SCHEME:' + name,
                            'data:image/svg+xml;base64,' + base64.b64encode(svg.encode()).decode())
    html = html.replace('/*ASDATA*/',
                        read(os.path.join(DATA, 'activeset_data.json')).strip())
    html = html.replace('/*REFDATA*/',
                        json.dumps(slim_reference(), separators=(',', ':')))
    html = html.replace('/*DATA*/',
                        json.dumps(json.loads(read(SITE)),
                                   separators=(',', ':')))
    html = inline_renders(html)
    html = html.replace('/*RELEASEDATA*/', json.dumps(json.loads(read(os.path.join(ROOT, 'pilot', 'joint_state_release', 'page_data.json'))), separators=(',', ':')))
    html = html.replace('/*DEPTHDATA*/', read(os.path.join(ROOT, 'web', 'depth_profiles.json')))
    html = html.replace('/*TOFDATA*/',
                        json.dumps(json.loads(read(TOF)), separators=(',', ':')))
    for token, path in PARTS:
        html = html.replace(token, read(path))
    for token, _ in PARTS + [('/*DATA*/', None), ('/*ASDATA*/', None),
                             ('/*REFDATA*/', None), ('/*TOFDATA*/', None),
                             ('/*RENDERDATA*/', None)]:
        if token in html:
            raise SystemExit('unsubstituted token left in the page: ' + token)
    # Links out (the repository, reference DOIs) are allowed; loading
    # anything from outside is not.
    if '<script src=' in html or 'href="http' in html.replace(
            'href="https://github.com', '').replace('href="https://doi.org/', ''):
        raise SystemExit('the page must not reach outside itself')
    return html


def build():
    html = render()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as fh:
        fh.write(html)
    return OUT, len(html)


if __name__ == '__main__':
    path, size = build()
    print('%s written (%d kB)' % (os.path.relpath(path, ROOT), size // 1024))
