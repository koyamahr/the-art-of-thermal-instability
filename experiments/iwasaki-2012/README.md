# Iwasaki & Inutsuka (2012): front-classification prototypes

Target: Figures 5 and 7 of *Structure of dynamical condensation fronts in the interstellar medium*, MNRAS 423, 3638–3645. The source inspected is arXiv:1204.5532v2 (14 May 2012); the PDF checksum and exact scope are in `source-audit.json`.

This is a **computed prototype on a restricted branch**, not a completed reproduction of all colored regions or auxiliary curves. The reference examples reproduce the nine classifications stated in Figures 2–4 and section 4.1.1. The map boundaries have not been digitized for a quantitative comparison.

## Run

Keep `front_model.py`, `reproduce.py`, `config.json`, `source-audit.json`, and `requirements.txt` together:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduce.py --output output
```

NumPy 2.0.2, Matplotlib 3.9.4, and SciPy 1.13.1 are pinned. This is a steady ODE calculation, requiring no 2D/3D hydrodynamic simulation. The validation report records the measured execution time for each environment; it is not a portable performance guarantee. Grid counts and tolerances can be changed in the configuration. In the repository, `publish.py` recomputes and synchronizes the public assets.

## Equations and units

The gas is pure atomic H with gamma=5/3. Use the source heating rate Gamma=2e-26 erg s^-1 per H atom, the cooling ratio of Eq. (7), and kappa=2500 sqrt(T) erg cm^-1 K^-1 s^-1. Define theta by dtheta=kappa dT and F=dtheta/dx. F is positive when temperature rises towards the warm side; the physical conductive heat flux is -F.

Mass and momentum conservation give j=rho*v and M=rho*v^2+P. The solver uses the **plus-density branch** of Eq. (11). Its temperature-form phase equation is

```text
B(T) = (R*j/2) * [(gamma+1)/(gamma-1) + M/sqrt(M^2-4R*j^2*T)]
d(F^2)/dT = 2B(T)*F + 2*kappa(T)*rho*L
```

The CNM boundary is its thermal equilibrium with F=0, approached along the outgoing saddle eigenvector. Matching independently integrated CNM and WNM branches at the unstable equilibrium yields j_ZP. The zero-flux saturation condition is Eq. (16), integral kappa*rho*L dT=0.

The numeric positive plotted flux is `-j / [m_H * (1 cm^-3) * (1 km s^-1)]`. Physical j has units g cm^-2 s^-1. The source axis abbreviates its normalization; the density unit is explicit here. The critical temperature is M^2/(4Rj^2), and the corresponding pressure is M/2.

For real eigenvalues at the target equilibrium, integrate log(F) and refine the remaining temperature distance. Capture gives F proportional to that distance, while crossing gives a finite F. Two consecutive distance refinements must agree. This follows the trajectory instead of treating a local node discriminant as a sufficient classification; the source's footnote 1 explicitly distinguishes these.

## What the prototype does not cover

- The source Figures 5/7 include regions for which this implementation's Eq. (11) branch has no outgoing CNM saddle. Such inputs are gray, not filled from a guess or copied from the paper.
- Numerically ambiguous cells stay white. Colors are computed only at the exported grid points; finite plotting cells are not precise boundary locations.
- Black denotes below/on the condensation ZP separatrix. Evaporating ZP solutions are outside this condensation-only calculation.
- Shock continuation, alternative density branches, multidimensional evolution, and the published j_crit overlay are absent.
- Section 4.1.2's choice of reference temperature for j_crit needs clarification, including the temperature interval and loss normalization. The direct trajectory classification does not require that overlay.
- Modern constants and the printed cooling fit give Pmin=1597.45, Pmax=5007.00 and Psat=2822.88 K cm^-3. The paper quotes Pmin=1596 in the text/1598 in Figure 1, Pmax=5012, and Psat=2823. These differences are retained.

The open review points describe implementation/interpretation questions, not an assertion of errors in the paper. A scientific review with the authors is the next step before calling the map a complete reproduction.

## Validation and output

Checks cover thermal landmarks, saturation quadrature, matched ZP flux, tighter solver tolerances, nine published classification examples and their refinement, mass/momentum conservation, critical-point pressure, and independent finite-difference derivatives. The validation report lists every tolerance and measurement. Its source/config hashes link it to the downloadable calculation.

Outputs include two classification PNG/SVG files; grid CSVs with class codes, physical fluxes, terminal temperatures, and decision reasons; boundary CSVs; an equilibrium PNG/SVG and CSV; and `validation.json`. The two grid CSVs are the underlying sampled results, including masked and unresolved cells.

Original implementation code is MIT. Generated data and figures have separate content terms; see https://koyamahr.github.io/rights/ . Original paper rights are unchanged.
