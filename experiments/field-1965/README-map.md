# Field (1965), Figure 1: domains of instability

Equation-based reconstruction of Figure 1 on printed p. 536, with definitions and boundaries checked visually against pp. 535–537 of the ADS scan. The range is -8 to 8 in both coordinates and gamma=5/3. Source provenance is in `../../sources/manifest.json`.

## Run

From the repository root, with the Figure 2 environment installed:

```sh
.venv/bin/python experiments/field-1965/reproduce-map.py --output experiments/field-1965/output
```

Standalone: keep `reproduce-map.py`, `config-map.json`, and `requirements.txt` together, install the requirements, then run `python reproduce-map.py --output output`. There are no local helper imports or PDF dependencies. The initial environment is Python 3.9.6, NumPy 2.0.2, and Matplotlib 3.9.4; complete package versions are recorded in `environment.txt`.

## Coordinates and cubic

Equation (17) defines `y=n/(k*c)`, `r=sigma_rho=k_rho/k`, and `t=sigma_T_prime=k_T/k+k/k_K`. Conduction is already incorporated in the horizontal coordinate. Notice that y differs from `x=n/(k_rho*c)` in Figures 2 and 3.

Equation (18):

```text
y^3 + t*y^2 + y + (t-r)/gamma = 0.
```

Define `A=t-r` and `B=(gamma-1)*t+r`. Away from the neutral lines the number of roots with **positive real part**, counted with multiplicity, is `(A<0) + 2*(B<0)`. This formula counts oscillatory growing roots as well as positive real roots.

- `r=t`: isobaric neutral line; there is a zero root.
- `r=-(gamma-1)*t`: isentropic neutral line; the cubic factors as `(y+t)*(y^2+1)`.

The lines are boundaries for the total count of growing roots. In the three-real-root regions they must not be described as boundaries of a currently oscillating wave: the acoustic roots there are non-oscillatory.

## Hatched regions are different boundaries

Field's equation (19) uses the Cardano quantity

```text
D = [(t-r)/gamma - t/3 + 2*t^3/27]^2/4 + [1-t^2/3]^3/27.
```

**D<0 gives three real distinct roots; D>0 gives one real root and a complex conjugate pair.** This sign convention is opposite to the usual cubic discriminant, which equals -108D. The two thorn-shaped regions are bounded by equation (20):

```text
r = (1-gamma/3)*t + (2*gamma/27)*[t^3 +/- (t^2-3)^(3/2)],
|t| >= sqrt(3).
```

The cusp points have `t=+/-sqrt(3)`, `r=(1-gamma/9)*t`, and the triple root `y=-t/3`.

## Connect to the growth-rate figures

Equation (28), with `q=k/k_rho>0`, gives `r=1/q` and `t=alpha/q+beta*q=alpha*r+beta/r`. The six conducting trajectories use alpha=1/2 and -2, and beta=0.01, 0.1, 1. Dotted beta=0 reference rays are also shown. Increasing q follows a trajectory toward smaller positive r; conduction shifts t to the right.

Mapping back to the growth-rate cubic gives `x=q*y`. The validation compares the entire root set along each path to the Figure 2/3 cubic and checks its intersections with the relevant neutral line. Thus the three exhibits use consistent normalizations.

## Validation and limitations

Checks compare analytical growing-root counts and D-based topology against numerical roots on a separate off-boundary grid, compare all roots against an independently assembled fluid matrix, evaluate scaled polynomial residuals, test the thorn boundaries and triple-root cusps, test neutral-line roots, and test path mapping and cutoff intersections. Any failure stops generation. The independent grid excludes points within 1e-6 of a neutral line or D=0.

Outputs: SVG and PNG figure, a 161x161 domain grid CSV, a trajectory CSV, and `validation-map.json` with measured errors, tolerances, environment, and source hashes. In the grid CSV, growing-root counts on neutral lines are blank; D-based topology uses an absolute 1e-10 tolerance for repeated roots. The plot shows only the requested window; paths outside it are clipped and omitted from the path CSV.

Root-count background colors and typography are modernized. The original circular guide for placing root-count labels is deliberately omitted; it is not a physical boundary. No digitized pointwise comparison with the original printed figure is claimed. Licenses remain unselected.
