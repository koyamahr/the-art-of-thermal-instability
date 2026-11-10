# Field (1965), Figure 2

An equation-based reproduction of the unstable condensation branch in Figure 2 of George B. Field, *Thermal Instability*, ApJ **142**, 531–567 (1965), DOI: [10.1086/148317](https://doi.org/10.1086/148317).

## Scope and source check

The published scan was downloaded from ADS; see `../../sources/manifest.json`. Equations (15)–(16) on printed p. 534, (27) on p. 537, (38)–(39) on p. 539, and Figure 2 on p. 540 were checked against the rendered pages. The configuration is alpha=1/2, gamma=5/3, beta=0, 0.01, 0.1, 1; horizontal range 0–7. These values match the labels and caption of the original figure.

This is a fresh calculation from the published equation. Original curve pixels have not been digitized, so there is no claimed pointwise error relative to the printed curves. Original typesetting and small scan artifacts are not reproduced.

## Derivation

Field defines a perturbation proportional to exp(n t + i k·r), so positive real n denotes growth. The heat-loss function L is losses minus gains **per unit mass per unit time**, excluding thermal conduction.

- `c = sqrt(gamma p0/rho0)`: adiabatic sound speed.
- `k_rho = mu (gamma-1) rho0 L_rho / (R c T0)`.
- `k_T = mu (gamma-1) L_T / (R c)`.
- `k_K = R c rho0 / [mu (gamma-1) K]`, with thermal conductivity K.
- `alpha = k_T/k_rho`, `beta = k_rho/k_K` (equation 27).
- `q = k/k_rho`, `u=q^2`, `x=n/(k_rho c)` (equation 38).

Here L_rho is at fixed T and L_T at fixed rho. For this figure k_rho>0. Beta=0 means the zero-conductivity limit.

Equation (39), regrouped as a cubic in x:

```text
x^3 + (alpha + beta*q^2)*x^2 + q^2*x
    + q^2*(alpha - 1 + beta*q^2)/gamma = 0.
```

For 0<alpha<1, beta>=0 and q>0 inside the unstable band, P(0)<0 and P'(x)>0 for x>=0. There is exactly one positive root. Bisection on `[0, (1-alpha)/gamma]` selects it without mode switching. At q=0 the limiting growth is zero. For beta>0 the neutral cutoff is `q_c=sqrt((1-alpha)/beta)`; beyond it the plotted branch is omitted and CSV growth values are blank. The full linear-system eigenvalues are still used to verify stability across the cutoff.

## Reproduce

From the repository root, using Python 3.9 or later:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r experiments/field-1965/requirements.txt
.venv/bin/python experiments/field-1965/reproduce.py --output experiments/field-1965/output
```

Alternatively download `reproduce.py`, `config.json`, and `requirements.txt` from the exhibit into a folder, install the requirements, and run `python reproduce.py --output output` there. No PDF, network access, random seed, or simulation inputs are needed after installing dependencies. The computation is a small CPU calculation. Output contains SVG, PNG, CSV, and a machine-readable validation report. `environment.txt` records the complete dependency versions of the initial successful run; `requirements.txt` pins direct dependencies.

Repository maintainers can run `.venv/bin/python experiments/field-1965/publish.py` from the repository root to recompute and synchronize the figures, data, validation report, and downloadable code into `docs/field-1965/`.

## Validation

The script fails if any check fails. It compares roots against an independently assembled linear system from equations (11)–(14), checks scaled polynomial residuals, tests the long-wave and short-wave analytical limits, and checks each conductive cutoff and the change from growth to damping across it. The tolerances and measured errors are stored in `validation.json`.

The linear system uses state `(delta rho/rho0, i delta v/c, delta p/p0)` and time `tau=k_rho*c*t`. With `s=alpha+beta*q^2`, its matrix is:

```text
[[0,     -q,       0],
 [0,      0, q/gamma],
 [s-1, -gamma*q,  -s]]
```

The limiting relations are `x/q -> sqrt((1-alpha)/(gamma*alpha))` as q→0 and, for beta=0, `x -> (1-alpha)/gamma = 0.3` as q→∞. Reported peaks are sampled maxima, not optimized estimates; beta=0 has no finite-wavenumber peak in this model.

## Publication status

Implemented and validated locally and in GitHub Actions. Original implementation code is MIT. Text, figures and generated data have separate content terms; see https://koyamahr.github.io/rights/ . Original paper rights are unchanged.

## Figure 3 is also available

See `README-wave.md` and `reproduce-wave.py` for the complex wave pair at alpha=-2. To rebuild only that exhibit, run `.venv/bin/python experiments/field-1965/publish.py --figure 3`. The default publish command rebuilds all three figures.

## Figure 1: the stability map

See `README-map.md` and `reproduce-map.py` for the root-count domains, discriminant boundaries, and conductive trajectories linking Figures 2 and 3. Run `.venv/bin/python experiments/field-1965/publish.py --figure 1` to rebuild the map only.
