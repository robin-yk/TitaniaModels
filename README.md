# TitaniaModels

Gas–solid equilibrium and oxygen-vacancy models for reduced titania.

Two calculations used in the accompanying manuscript: gas–solid equilibrium for
the C–H–O–N / Ti–O system, and the isolated surface-vacancy population that sets
the CO formation rate across the reduction series. The page also shows a pilot
vacancy-distribution model that the manuscript does not use.

The interactive version is at
[robin-yk.github.io/TitaniaModels](https://robin-yk.github.io/TitaniaModels/).

## What it computes

### Gas–solid equilibrium

Given a gas charge and a temperature, Gibbs-energy minimisation returns the
stable Ti–O phase assemblage, the equilibrium gas composition, the CO2
conversion, and the margin to the nearest reduced phase. Condensed phases are
carried as an active set, so a phase that is not stable has exactly zero moles
rather than a numerical floor, and its distance from stability is reported as a
KKT reduced cost instead. The feed at which rutile stops surviving comes out in
closed form.

Gas data are NIST-JANAF; the titanium oxides through Ti2O3 are the Waldner and
Eriksson CALPHAD assessment.

The equilibrium workspace also carries a thermodynamic operating atlas.  It
first solves fresh-feed RWGS without a solid to show the gas-phase conversion
ceiling across temperature and H2:CO2 ratio.  It then overlays the TiO2
reduction boundary and repeats a ratio sweep with the finite oxide charge.  The
two calculations are kept separate so gas reaction equilibrium is not confused
with oxygen supplied by the solid.

### Pilot: time-dependent vacancy distribution

The distribution tab shows the R600 reduction and reoxidation calculation in
`pilot/joint_state_release/`. It follows vacancies and Ti³⁺ from the surface
to the centre of a 900 nm particle at 600 °C. The calculation starts from a
defect-free particle, uses the measured reduction history, and then follows
600 s of reoxidation. The selected mobility case uses assumed barriers of
1.6 eV in the outer 30 nm and 2.2 eV below that region. The page labels these
inputs and displays stored results. Source code, full results, and provenance
are included in the release folder.

The manuscript does not use this pilot. Its Supplementary Note 2a divides every
sample's initial CO rate by one fixed reactive-site count, 2.31 umol/g, and its
Note 2b profile comes from a continuum model with local electroneutrality that
is no longer in this repository. The pilot model in
`pilot/titania-super-multiscale/` replaces that number with a range. At
equilibrium at 873.15 K with the measured inventory fixed, it places each
sample's vacancies over every site a vacancy can occupy: explicit (110)
trilayers with bridging, in-plane and sub-bridging oxygen, (1×2)
reconstruction cells, Monte Carlo aggregate boxes in the bulk, and shell
electrostatics with only the whole particle neutral, over an equal-mass
900–1600 nm size mixture. Values without a literature source are the model's
parameters. The model is solved at 720 parameter points under five
reactive-site definitions, and the rate divided by the reactive-site count
over the countable cases (at least 0.01 umol/g and 1% of the bridging
capacity) is the apparent TOF range; the other cases are kept and reported
apart. The run takes about 40 minutes on four cores.
`scripts/export_tof_range.py` packs it into `paper_outputs/tof_range.json`,
and the historical export remains available for reproduction.

### Vacancy population model

Titration gives the total oxygen removed per sample; it rises eightfold across
the reduction series while the initial CO rate peaks after 600 °C. The smallest
model that produces that shape counts only isolated vacancies in a surface
active region of capacity nₛ, generated uniformly during reduction and lost by
a thermally activated second-order association. Writing q = n_iso + n_assoc the
loss term cancels from the sum, so q and n_below are closed form and one scalar
equation is integrated; the mass balance is an identity. Three parameters
(nₛ, E_loss, ν_eff) are fitted to five rates, and the residual definition is
reported with every fit because the rates span 86-fold.

## Reproducing the manuscript numbers

```bash
python3 scripts/reproduce_paper.py    # -> paper_outputs/*.csv and site_data.json
python3 scripts/build_site.py         # -> docs/index.html
python3 -m pytest tests/ -q
```

`paper_outputs/` is committed, and a gate fails if regenerating it changes a
byte. The same gate holds `docs/index.html`; its MathML depends on the
latex2mathml version, which `requirements.txt` pins.

## Verification

Three implementations must agree before a number is quoted: the Python package,
an mpmath oracle at 80 significant digits in `scripts/oracle_tio.py` that shares
only the raw coefficient tables and never imports the solver, and the browser
mirror. The page calculates results for the entered composition and parameters.
Parity tests compare each browser engine with its Python implementation:
`web/activeset.js` to 4 ulp (the difference between glibc's `exp` and
V8's), `web/population.js` to a relative 1e-6, the difference between two
integrators of the same equation.

The fresh-feed gas-only RWGS relation has its own smaller 80-digit oracle in
`scripts/oracle_rwgs.py`.  Its temperature-ratio map is evaluated in the browser
from the coefficient tables.

The page loads manuscript reference values from committed JSON and keeps them
separate from calculations using user inputs.

The historical TOF range has no browser engine. Gates hold
`paper_outputs/tof_range.json` to the pilot outputs byte for byte.
The current distribution tab reads the joint-state release data.

`scripts/check_page.py` loads the built page in a browser and compares its
displayed numbers with fresh calculations.

GitHub Actions (`.github/workflows/gates.yml`) runs the test suite, the pilot
tests and the browser check on every push and pull request, and fails on any
skipped or xfailed test.

## Layout

```
solidgas/       activeset · vacancy_population · species · shomate · waldner
analysis/       dielectric_correlation — an empirical regression, not thermodynamics
scripts/        reproduce_paper · build_site · oracle_tio · oracle_rwgs · check_page
data/           thermodynamic tables, the reduction series, the 80-digit reference
paper_outputs/  every computed number in the manuscript, as committed CSV
web/            the browser mirrors, a drawing kit, figure modules, one page
docs/           the built page and the method notes
pilot/          joint_state_release; historical titania-super-multiscale TOF model
```

See [docs/equilibrium-method.md](docs/equilibrium-method.md) and
[docs/verification.md](docs/verification.md).
