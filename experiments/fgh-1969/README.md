# FGH/GHF 1969: model ingredients

This is an **ingredient diagnostic**, not a reproduction of FGH's equilibrium Figures 1 and 2. It evaluates published heating expressions and cooling-transition fits, supplying inputs for the separate reduced cold-branch experiment described below.

## Sources

- Field, Goldsmith & Habing (1969), *Cosmic-Ray Heating of the Interstellar Gas*, ApJ 155, L149–L154: the target equilibrium figures and assumptions.
- Goldsmith, Habing & Field (1969), *Thermal Properties of Interstellar Gas Heated by Cosmic Rays*, ApJ 158, 173–184: detailed companion paper, identifying the Letter as Paper I. Implemented: Table 1 and equations (16), (23), (25), (26), on pp. 178–181.
- Field et al. (1968), ApJ 151, 953: earlier cooling data; downloaded as a supporting source but not fully audited.

Table 1 and the implemented equations were visually checked against rendered scans. Bibliographic source URLs and PDF hashes are in `../../sources/manifest.json`; `source-audit.json` records the implemented pieces and unresolved inputs.

## Run

Use the same pinned NumPy/Matplotlib environment as the Field exhibits, or install this folder's requirements. From the repository root:

```sh
.venv/bin/python experiments/fgh-1969/reproduce.py --output experiments/fgh-1969/output
.venv/bin/python experiments/fgh-1969/publish.py
```

For a standalone download keep `reproduce.py`, `config.json`, `source-audit.json`, and `requirements.txt` together. Install the requirements and run `python reproduce.py --output output`. No original PDF or helper module is needed. Outputs include two figures (PNG/SVG), two CSVs, and `validation.json` with the actual runtime and source hashes.

## Heating definitions

Let `x=n_e/n(H I)`, `n_target=n(H I)+n(He I)+n(He II)`, and `Q=(1+130*x)/(1+57*x)`.

- Equation (26): `e_h=9*Q` eV per primary liberated electron.
- Equation (25): `f=0.13*Q`, a fraction of the cosmic-ray energy loss to **ionization plus excitation**, not a fraction with the same denominator as e_h.
- Equation (23): `Gamma_ion/(zeta*n_target)=1.4e-11*Q` erg.
- Direct electron heating, stated on pp. 177 and 179, is `Gamma_el=4.6e-10*n_e*zeta` erg cm^-3 s^-1. It is not added to the plot without specifying the particle populations.

These are historical rounded formulas. Equation (26) gives approximately 9, 13, and 20 eV for x=1e-4, 1e-2, and 1, as stated in the paper; its printed fit slightly exceeds 20 eV near x=1. The code retains that result and does not extrapolate the expression as a physical energy bound. Likewise, the coefficient 0.13 in equation (25) is not identical to 1.4/11 from the rounded coefficients of equations (23)–(24).

## Cooling definitions

Equations (15)–(16):

```text
r_k(T)=R*T^m*exp(-T_star/T)
Lambda_2=sum_k n_k*(n_e*r_ek + n_H*r_Hk)
```

The r values are energy-loss coefficients in erg cm^3 s^-1, not bare collisional excitation probabilities. The tabulated R has additional K^-m units. Five electron fits and six H-atom fits are available. Blank entries stay null and cannot be evaluated as a measured zero.

The plot evaluates the historical fits from 10 to 20,000 K and clips the displayed y range below 1e-29. This chosen evaluation range does not establish the accuracy of every fit throughout it. CSV values store log10(r), preserving coefficients that would underflow in ordinary floating-point exponentiation. The curves are coefficients per collider/target pair: densities, actual level populations, and historical abundances are needed before interpreting relative contributions to total cooling.

## What remains before equilibrium figures

The short Letter does not contain all rate details. The detailed paper requires coupled balances for temperature, electron density, superthermal electrons, and H/He ionization, plus historical atomic data and abundances.

There is also an unresolved input discrepancy: the Letter's Figure 1 caption gives Fe depletion factors 4 and 10 for cases a and b. The detailed paper's section IV gives remaining Fe fractions 0.3 and 0.01, corresponding to depletion factors about 3.33 and 100. The source audit preserves both statements. They must not be silently combined into a supposedly exact reproduction.

Validation checks the quoted heating anchors, limits, monotonicity, efficiency bounds, analytical cooling slopes, missing-channel handling, and the C II collider-contribution ratio. Passing these checks verifies the implemented ingredients, not the full FGH equilibrium solution. Licenses remain unselected; original paper rights are separate.


## Reduced cold-branch experiment

Run `python cold_branch.py --output output` with `cold_branch.py`, `reproduce.py`, `config.json`, `config-cold.json`, and `requirements.txt` in one folder. `publish.py` regenerates and publishes both experiments. Outputs: `cold-branch.svg`, `.png`, `.csv`, and `validation-cold.json`.

This is **not a reproduction of FGH Figures 1 or 2**. It solves primary-only H/He ionization and charge balance at each trial temperature, then scans for thermal-balance roots within 10–2,000 K. Root brackets are refined by bisection, and a doubled scan checks the numerical roots. A sign-change scan can miss tangential roots. Missing segments mean no bracket was found in this window, not that no physical equilibrium exists.

Hjellming (1966), p. 435, equations (56)–(59), gives total and ground-state recombination number (`beta`) and energy (`gamma`) coefficients. We subtract the ground-state fits for excited-state recombination and use `h*nu0*(gamma-gamma1)` for kinetic-energy cooling (definition on p. 425). The exponent of lambda in the last term of gamma is **−1/3**. He I is approximated by H, and hydrogenic Z=2 scaling supplies He II recombination. The paper's cited accuracy for the He approximation concerns H II region temperatures; it does not establish cold-gas accuracy.

The baseline abundances come from Field et al. (1968), Table 1, p. 958, citing Allen (1963): per million H nuclei, C=300, O=670, Si=30, Fe=8. We use the Letter's depletion divisors explicitly, retaining the different companion-paper values in the audit. Dashed curves use H/He charge only; solid curves add one electron per C, Si, Fe nucleus. Other trace donors are absent. All transitions use the elemental target density in a low-density approximation, without a statistical-equilibrium level calculation. Blank collider entries are omitted channels, not measured zeros.

Heating is GHF equation (23) times `zeta*n_target`, plus direct electron heating. Recombination energy follows Hjellming; free-free cooling follows GHF equation (18). Secondary ionization, thermal H/He ionization and excitation, molecules, grains, and stellar heating are absent. The pressure proxy `(n0+ne)*T` omits trace nuclei. No stability assignment, warm branch, or exact agreement with the original Letter is claimed. The config and validation report carry these assumptions.

The previous ingredient CSV mislabeled its normalized heating denominator as `n0`; it is corrected to `n_target`. Its normalized numerical values are unchanged.


## Cold-to-warm equilibrium prototype

Run `python two_phase.py` beside `cold_branch.py`, `reproduce.py`, `config.json`,
`config-cold.json`, and `config-two-phase.json`, with the pinned requirements.
It writes `output/two-phase.{png,svg,csv}` and `validation-two-phase.json`.
No SciPy or source PDF is needed. Parameters and approximations are recorded in
config-two-phase.json; changing code constants requires a corresponding input audit.

The calculation solves charge neutrality, H/He ion populations, superthermal
production/removal, secondary hydrogen ionization, and energy balance. It uses
GHF Eq13 directly rather than combining its rounded normalization with Eq23.
The 20-eV excitation point is inferred from Eqs12-13 and the explicitly printed
H ionization triangle. A linear excitation ramp from 10.2 to20eV is assumed;
its shape and normalization have not been validated against the original
Fite et al. data or Goldsmith thesis. Zero and continued linear tails beyond20eV
are sensitivity alternatives, not physical uncertainty bounds.

Temperature parameterization retains folds in T(n). A sign-change density scan
and doubled scan are compared; roots outside0.01–1000cm^-3 or10–12000K are
not extrapolated. The report lists sampled temperatures without a bracket.
Pressure turning points are grid brackets for this prototype, not historical
measurements, conductive saturation pressures, or a full stability analysis.

The baseline uses Hjellming kinetic recombination energy. The separate2/3 test
records the unresolved GHF Eq17 convention. Thermal He excitation/ionization,
metal level populations, and other limitations remain explicit. No modern KI or
Wolfire fit substitutes for the historical H cooling here. Passing48 numerical
checks is not a claim to reproduce the original FGH figure quantitatively.


## Original-figure comparison and stopping note

Historical coefficient recovery is **deferred**. The user explicitly permits
moving on when exact reproduction has no clear prospect. No coefficients have
been fitted to the original curves. The existing prototype remains unchanged.

Run `python compare.py` with all existing prototype dependencies, plus
`digitized-fgh.csv` and `config-comparison.json`, in the same directory.
The script writes `output/comparison.{png,svg,csv}` and
`output/validation-comparison.json`. A source PDF is not needed to rerun the
comparison. The JSON identifies its SHA256, page indices, render scale,
crop origins, axis tick coordinates, and reading method.

The 96 sparse readings were selected on isolated curves and checked visually
against the source. Pixel envelopes include line width, a practical3-pixel
reading envelope and axis-fit residuals; they are not statistical confidence
limits. No independent second digitization has been performed. Source scans
and the private overlay used for checking are not distributed.

Figure1 compares density at each reference temperature. Figure2 compares
pressure/electron density at the reference density, retaining every temperature
root; ambiguous or absent roots are recorded rather than choosing the closest
branch. This sample has no unscored points. The metrics differ and must not be
pooled into a single reproduction score. Refining the scans preserves the roots.

The refined model warm pressure maximum is1516.695Kcm^-3 at n0=0.165636cm^-3,
T=8417.37K. The text quotes1800Kcm^-3,0.21cm^-3,7800K. The pressure discrepancy
is-15.74%, density-21.13%, temperature+7.91%. At the quoted n,T point, model
cooling/heating=1.17138. The shifted rapid transition makes some same-density
pressure discrepancies much larger: the largest sampled difference is1.10dex.
It would be misleading to describe the entire figure as agreeing within16%.

The depleted a/b curves in the original are cold-side segments; the caption
assumes transition to normal abundances near n~1cm^-3. The comparison restricts
the plotted model segments accordingly. The separate full-domain prototype
shows constant-abundance model extensions, not unseen original curves.

### Coefficient audit

- GHF Table1 (p179), Eq13 and16 (p178), and Field1968 Table1 (p958) were visually
  rechecked. No transcription error was found in the currently used metal fits
  or abundance baseline. This does not establish that every approximation
  matches the unpublished original implementation.
- Field1968 Eqs17–18 discuss Fe cascades. The GHF554/961K exponential terms
  should not automatically be reinterpreted as independently populated lower
  levels and suppressed by additional level fractions. The low-density total
  cooling sum is retained; finite-density level populations remain unvalidated.
- At the printed n,T anchor, SiII contributes31.4% of model cooling and
  OI22500 contributes18.0%; all term budgets are in the report.
- Eq13+20 reproduces rounded Eq23 heating within3%. This rounding difference
  is not a resolved explanation of the pressure discrepancy.
- GHF Eq17 prints a2/3 factor with its energy-weighted recombination definition.
  The baseline retains the Hjellming kinetic-energy convention; a separate
  printed-factor probe moves the warm pressure maximum by less than1%.
- Halving/doubling thermal Ly-alpha cooling moves the warm maximum by about
  +1.6%/-1.9%. This is a term-sensitivity probe, not a recovered historical
  cross section or a bound on all H/He atomic-data uncertainties.

The original Ly-alpha fit/Goldsmith thesis, thermal He processes, recombination
convention, and exact historical implementation remain unresolved. The printed
FGH and GHF depleted cases also differ and are not merged. Resume only if a
specific missing fit, original code, or other discriminating evidence becomes
available; continued unconstrained coefficient adjustment is not planned.

The80 comparison checks test calibration, roots, consistency, and budgets. They
do not turn the measured figure discrepancies into validation passes.

## Supplementary bounded ZP search

After two_phase.py, run `python saturation.py`. See config-saturation.json.
This adds an assumed Parker conductivity and instantaneous ionization closure.
Normal abundances: no sign bracket at the sampled pressures inside the retained
10--12000 K domain; P_sat is null. Depleted a/b warm extensions have conditional
roots, stored with all three temperatures, H+He nuclei density n0, hydrogen nuclei
density 0.9 n0 and electron density. These are not historical FGH saturation values.
The finite search does not exclude an unsampled or out-of-domain root.
Historical coefficient recovery remains deferred.
