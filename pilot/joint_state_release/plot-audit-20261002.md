# Distribution plot audit, 2 October 2026

## Defect and correction

Oxygen-only atomic planes have zero titanium capacity. Their Ti³⁺ fraction is undefined. The source exporter uses zero as a placeholder. The plot previously changed this placeholder to 10⁻¹⁰% for the logarithmic axis.

The plot now excludes points with zero site capacity. It also excludes non-positive values from logarithmic axes. Valid atomic-plane points keep their original indices. A missing point breaks a bulk line. Source calculations and saved numerical results are unchanged.

## Checks

- Checked 3,016 saved vacancy and Ti³⁺ fractions across all saved times. Each defined fraction equals 100 × amount / site capacity.
- Checked non-negative amounts and site capacities, and amounts no greater than capacity.
- Checked total electron amount = 2 × total vacancy amount at every saved time, within 10⁻⁶ µmol g⁻¹.
- Found 104 undefined Ti³⁺ entries on oxygen-only planes; all have zero electron amount.
- Checked all 42 reduction/reoxidation time pairs in the browser. SVG coordinates are finite. Atomic-plane point counts match the positive-capacity data.
- The existing browser gate passed, including gas equilibrium parity, population mass balance, four release figures, saved endpoints and equation rendering.
- Paper outputs regenerated without changes to tracked files.
- Checked exact source parity for 1,964 exported fields.
- Full test suite: 239 passed. One overflow warning occurred in the gas-equilibrium active-set test; that test passed.

These checks test numerical consistency and plotting. Transport barriers remain the stated scenario inputs.
