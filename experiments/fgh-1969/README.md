# FGH/GHF 1969: model ingredients

This is an **ingredient diagnostic**, not a reproduction of FGH's equilibrium Figures 1 and 2. It evaluates published heating expressions and cooling-transition fits, supplying explicit inputs for a future coupled equilibrium solver.

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

Let `x=n_e/n(H I)`, `n_0=sum_j n_j`, and `Q=(1+130*x)/(1+57*x)`.

- Equation (26): `e_h=9*Q` eV per primary liberated electron.
- Equation (25): `f=0.13*Q`, a fraction of the cosmic-ray energy loss to **ionization plus excitation**, not a fraction with the same denominator as e_h.
- Equation (23): `Gamma_ion/(zeta*n_0)=1.4e-11*Q` erg.
- Direct electron heating, stated on pp. 177 and 179, is `Gamma_el=4.6e-10*n_e*zeta` erg cm^-3 s^-1. It is not added to the plot without specifying `n_e/n_0`.

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
