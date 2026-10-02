# R600 joint-state transport release

Model ID: `model-tio2-surface-bulk-joint-state`, version `0.1.0`.

This release preserves the R600 reduction and reoxidation calculation from the wiki depth-barrier study dated 2026-09-29. The source files and full result arrays are copied without numerical changes. `provenance.json` gives the SHA-256 hash of each source file.

## Selected calculation

The particle diameter is 900 nm. The temperature is 873.15 K. The outer 30 nm has a transport barrier of 1.6 eV. The interior barrier is 2.2 eV. The attempt frequency is 10¹³ s⁻¹ and the hop distance is 0.3 nm. The model starts without defects, applies the measured reduction rate for 1800 s, then applies the measured reoxidation rate for 600 s. The measured rates are boundary inputs.

This is a conditional transport test for R600. The barrier profile was selected with both measured histories. It does not predict either measured rate. It does not establish a unique physical barrier profile.

The reduction integral is 95.7192445 µmol O g⁻¹. The final bridging-site vacancy fraction is 39.95657%. The depth that contains 90% of the vacancy inventory is 37.205 nm. The imposed reoxidation removes 81.525 µmol O g⁻¹. The final bridging-site vacancy fraction is 3.19214%.

## Files

- `source/`: unchanged numerical source and required input files.
- `results/`: complete reduction, reoxidation, and case-summary JSON files.
- `page_data.json`: display export. It retains every scalar checkpoint. It retains all six reduction profiles and seven reoxidation profiles. Each display profile has at most 120 depth points, including the first 12 atomic planes. Full profiles remain in `results/`.
- `depth-dependent-mobility.svg` and `.png`: original study figure.
- `export_from_wiki.py`: copies the release from the wiki.
- `verify_release.py`: checks source hashes and display values against the full results.

## Reproduce

Use Python with NumPy and SciPy installed. From this directory:

```sh
python source/titania-dielectric-redox/models/time-dependent-redox-transport/studies/depth-barrier-test/2026-09-29/run_case.py --outer 1.6 --inner 2.2 --depth 30 --bulk 320 --reduction-step 2.5 --reox-step 0.05
python verify_release.py
```

The calculation writes new results beside `run_case.py` and refuses to replace existing outputs. The imported result has a historical reduction filename that says `n32`; its stored grid has 320 bulk cells. Use the JSON metadata for the actual grid size.

## Source paths

The original scripts use historical wiki paths. This source tree preserves those paths so the scripts can run unchanged. The wiki now stores `time-dependent-redox-transport` under `models/vacancy-transport/`, and `atomic-layer-and-bulk-defect-distribution` under `models/vacancy-distribution/`. The historical study names `surface-bulk-connection` and `depth-barrier-test` map to `joint-state` and `depth-barrier`. Canonical-path dependencies are also copied because one module uses the new paths. The source root keeps the name `titania-dielectric-redox`, which the scripts use to locate their inputs.

Raw results preserve their historical source metadata. `provenance.json` records the files exported with this release. It does not replace the historical metadata.

## Numerical checks

The saved study includes mass balance, electron balance, Jacobian, restart, positivity, and free-energy checks. Source loading was tested through `run_case.py --help`. The export check verifies byte-identical source hashes and exact scalar checkpoint values. The full calculation was not repeated during this import.

## Publication figure export

The four page panels export at 5 × 5 inches. Each inner axis is 262.8 × 262.8 points (3.65 × 3.65 inches). Helvetica is the first font, with Arial as fallback. Axis labels are 18 points, tick labels 15 points, and legends 14 points. Axis and curve widths are 2 points. The measured inventory uses a 12-point open circle. Each panel has an editable SVG export and a 3000 × 3000 pixel PNG export (600 pixels per inch at 5 inches). Browser checks confirmed the SVG dimensions and square axes. All four default panels were visually checked; depth-profile legends are below the curves.
