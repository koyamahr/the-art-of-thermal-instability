# Takayanagi & Nishimura (1960): source rates and a finite level ladder

Source: *Cooling of the Interstellar Clouds in Region of Neutral Hydrogen*, PASJ 12, 77–105, DOI [10.1093/pasj/12.1.77](https://doi.org/10.1093/pasj/12.1.77). Table 1 and equation (7), p. 82, and equation (9), photon energy, and Table 2, p. 83, were visually checked against the ADS scan. The PDF checksum is in `source-audit.json`. Source PDFs are not included in the publication.

This exhibit transcribes 34 **upward** excitation coefficients for H-atom collisions with para-H2, in cm^3 s^-1. They are not de-excitation coefficients or energy-loss coefficients. The temperature samples span 20–5,000 K; different transitions have different tabulated coverage. This is coverage of the source table, not a modern claim of physical accuracy over that range. No interpolation or extrapolation is performed.

## Run

Keep `reproduce.py`, `config.json`, `source-audit.json`, and `requirements.txt` together:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduce.py --output output
```

In the repository, `python experiments/takayanagi-nishimura-1960/publish.py` regenerates and synchronizes the public assets. Only NumPy and Matplotlib are required. The calculation uses at most five levels and needs no hydrodynamic simulation, external solver, or original PDF at runtime.

## Equations and normalization

For each adjacent even-J pair, Table 1 supplies `q_up`. With J the **upper** level, the historical expressions are:

```text
q_down = q_up * (2J-3)/(2J+1) * exp[170.68*(2J-1)/T]  cm^3 s^-1
A_down = 7.52e-13 * J*(J-1)*(2J-1)^4/(2J+1)          s^-1
DeltaE = 2.3561e-14 * (2J-1)                         erg
p_J/p_(J-2) = n(H)*q_up / [n(H)*q_down + A_down]
L_per_molecule = sum_J p_J * A_down * DeltaE          erg s^-1
```

Normalize with `sum p_J=1`. Multiplying `L_per_molecule` by the para-H2 number density gives a volumetric loss rate. `n(H)` is the H-atom collider density, not the total H-nuclei density. The rounded temperature and photon-energy prefactors are preserved independently as printed.

The model assumes steady, optically thin populations, H-atom collisions only, spontaneous downward emission, and only para-H2. At each temperature, close the ladder at the highest connected even-J level supplied by Table 1: Jmax=2 below 300 K, 4 at 300 K, 6 at 500 K, and 8 at 1,000–5,000 K. Missing entries remain null in the source configuration. Closing the finite model does not establish zero physical excitation to unlisted levels. It omits those channels as a declared approximation. No chemical or thermal balance of the gas, ortho-H2, line trapping, or other colliders are modeled.

## External comparison versus numerical validation

`table2-comparison.csv` compares 60 printed cooling values at four collider densities. The maximum difference is approximately -21.72% at T=200 K and n(H)=100 cm^-3: the truncated model gives 4.657765e-25 erg/s per molecule, versus the printed 5.95e-25. The cause is not established. Missing higher-level information, source precision, or a source error must be investigated before assigning a correction. Inputs are not fitted to force agreement. A source row appears as 155 K in Table 2 while Table 1 gives 150 K; it is excluded pending clarification.

The 12 numerical checks include an independent constrained rate-matrix solution, stationary fluxes, population normalization and positivity, equality of collisional energy loss and radiative emission, the dilute limit, and the dense Boltzmann limit of the **same finite ladder**. Table differences are recorded separately and are not labeled as passed checks. The source's p. 83 discussion warns about its potential approximation and high-temperature accuracy; no independent modern rate-data validation is claimed.

## Files

- `rotational-rates.csv`: 34 upward coefficients plus derived downward coefficients, radiative rates, and photon energies.
- `cooling.csv`: 64 finite-ladder solutions, Jmax, and level populations; unmodeled levels have blank cells.
- `table2-comparison.csv`: 60 source/model comparisons and fractional differences.
- `validation.json`: numerical checks, runtime, hashes, limitations, and separate external-comparison statistics.
- Two PNG/SVG figures: rates and cooling, plus discrepancies from the printed cooling table.

Original implementation code is MIT. Generated data and figures have separate content terms; see https://koyamarium.github.io/rights/ . Original paper rights are unchanged.
