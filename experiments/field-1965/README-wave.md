# Field (1965), Figure 3: wave mode

This calculation reproduces the growing oscillatory wave branch for **alpha=-2**, **gamma=5/3**, and **beta=0, 0.01, 0.1, 1**. Source: Field, *Thermal Instability*, ApJ 142, 531–567 (1965), Figure 3 on printed p. 540. The caption, curve labels, and axes were checked against the published scan. Equations (15), (27), (30), (38)–(40) supply the dispersion relation, normalization, and analytical checks.

## Run

Use the same environment as Figure 2 (Python 3.9.6, NumPy 2.0.2, Matplotlib 3.9.4 in the recorded run). From the repository root:

```sh
.venv/bin/python experiments/field-1965/reproduce-wave.py --output experiments/field-1965/output
```

For a fresh environment, use the setup instructions in `README.md`. For a standalone download, keep `reproduce-wave.py`, `reproduce.py`, `config-wave.json`, and `requirements.txt` together, install the requirements, and run `python reproduce-wave.py --output output`. The shared `reproduce.py` supplies the original cubic and the independent fluid matrix; its Figure 2 main function is not executed by the wave script. No PDF or network access is needed after installation.

## Select the wave, not the fastest root

Write `q=k/k_rho`, `x=n/(k_rho*c)`, and `s=alpha+beta*q^2`. The cubic is

```text
x^3 + s*x^2 + q^2*x + q^2*(s-1)/gamma = 0.
```

The script selects the root with **positive imaginary part**. Its conjugate has the same growth rate. For the four Figure 3 curves a unique complex pair exists throughout the sampled interval. If this is not true, the script stops instead of silently selecting another mode. At q=0 the pair has the limiting value x=0, while the condensation root is 2. Selecting the largest real part would therefore produce the wrong figure.

The regrouped equation (40) independently checks the real part xi:

```text
xi^3 + s*xi^2 + (s^2+q^2)*xi/4
     + q^2*((1-1/gamma)*s + 1/gamma)/8 = 0.
```

## Physical checks

For small q:

```text
Re(x)/q^2 -> -[1+(gamma-1)*alpha]/(2*gamma*alpha^2) = 0.025
Im(x)/q   -> sqrt[(alpha-1)/(gamma*alpha)] = sqrt(0.9)
```

The leading real-part coefficient is independent of beta. For beta=0 at large q, `Re(x)->0.1` and `Im(x)/q->1`.

From equation (26), the wave cutoff is

```text
q_c = sqrt[(-alpha-1/(gamma-1))/beta] = sqrt(0.5/beta).
```

At the cutoff `s=-1/(gamma-1)=-1.5`, the cubic factors as `(x+s)*(x^2+q^2)`. The wave pair is `x=+/-i*q`, but the third root is `x=1.5`: condensation is still unstable. Damping of the wave alone must not be described as stability of the entire gas.

Validation compares the selected roots and their conjugates with the fluid matrix, checks both polynomial residuals, checks long-wave and short-wave limits, checks neutral growth and frequency at each cutoff, checks growth changes to damping, and explicitly guards against substituting the fastest condensation root. Measured errors and tolerances are recorded in `validation-wave.json`. Any failure stops figure generation.

## Outputs and limitations

- `figure-3.svg` and `.png`: positive wave growth and neutral endpoints, q=0–7.
- `figure-3.csv`: full signed wave growth and positive frequency, including damped points omitted from the plot.
- `validation-wave.json`: parameters, checks, sampled peaks, environment, and hashes of both source scripts and configuration.

Overall curve shapes are visually checked. No digitized pointwise comparison with the original printed curves is claimed. Peaks are sampled, and the beta=0 maximum over the plotted range is an endpoint rather than a physical finite-wavenumber maximum. Original implementation code is MIT; content and data have separate terms, as for Figure 2.
