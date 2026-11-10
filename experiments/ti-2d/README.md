# Two-dimensional thermal condensation — cover illustration

This is a real, freely evolving compressible-flow calculation, not a historical-paper reproduction and not sustained turbulence. Heating/cooling uses the corrected KI-type analytic fit. Thermal conductivity is deliberately enlarged and made cold-enhanced to broaden interfaces on a modest grid. Do not use this movie to infer physical cloud sizes, growth times, or the Parker-conductivity saturation pressure.

## Model and units

Periodic 4.8 x 3.0 pc, 256 x 160 cells, gamma=5/3, m=1.37 m_H, n=rho/m, P=n k_B T. Mean n=2.8 cm^-3. Low-wave-number Fourier perturbations have density RMS 10%, seed 1965, zero velocity and initially uniform equilibrium pressure. The same continuous Fourier pattern initializes the 128 x 80 comparison.

Gamma=2e-26 erg s^-1 per particle; net volumetric cooling Q=n^2 Lambda(T)-n Gamma. Lambda/Gamma=1e7 exp[-118400/(T+1000)] + 0.014 sqrt(T) exp(-92/T), in cm^3. The first coefficient uses the corrected KI fit; see the archive's KI2002 source audit. Here n is explicitly rho/(1.37 m_H), with no separate chemistry or ionization solver.

kappa(T)=3e6 + 3e8 (100 K/T)^2 erg cm^-1 s^-1 K^-1. This illustrative transport law is **not** the Parker law or the original KI conductivity. Its own ZP pressure is computed independently and stored in diagnostics.json; it must not replace the existing KI/Parker benchmark.

## Numerical method

Conservative finite volumes; primitive-variable monotonized-central slopes; HLLC fluxes with an HLLE fallback if a star state is invalid; SSPRK2. Centered harmonic-face conduction and radiative source terms are advanced in the same RK stages. Cooling is log-tabulated at 16385 points over 2–30000 K; out-of-domain states cause a step retry, not extrapolation. Positivity failures reduce the global timestep. There are no density or temperature floors, clipping, imposed condensations, forcing or mass sinks.

Steps obey hyperbolic CFL 0.35, diffusion factor 0.18 and radiative fractional internal-energy change 0.12. Binary state files contain little-endian float64 U=(n,n vx,n vy,E), with x fastest, then y, and component last. n is in cm^-3, v in pc/Myr, E and pressure in m (pc/Myr)^2 units times cm^-3. The solver constants gamma and mass are fixed in source. Diagnostics use actual physical times. Frames are not reversed to fabricate a cycle.

## Reproduce

Requires a C++17 compiler, Python 3.9+ and requirements.txt. From this directory:

```sh
python -m pip install -r requirements.txt
python validate.py
python -c 'import json; from model import *; c=json.loads((HERE/"config.json").read_text()); b=compile_solver(Path("/tmp/ti2d-build")); run(b,initial(c),c,HERE/"output/raw",frames=c["snapshots"]); run(b,initial(c,128,80),c,"/tmp/ti2d-coarse",frames=24)'
python render.py
```

The renderer uses bundled ffmpeg to encode actual density snapshots. On systems lacking Helvetica it uses Pillow's default font. Exact video bytes may differ between ffmpeg/platform versions; physical states and diagnostics are the scientific outputs. No executable or raw frame sequence is distributed. Regeneration needs a few hundred MB of disk space. The top-level publish.py is a repository maintenance tool, not required to reproduce the downloads.

## Checks and limits

validate.py checks an analytically advected contact, continuum heat-equation decay, cooling against independent adaptive ODE integration, and a compressible linear TI eigenmode derived from an independent 3x3 Fourier matrix (density, velocity and pressure). It also checks pressure and energy conservation. `--check-only` runs benchmarks without replacing the archived full-run report.

history.csv stores global integrals of solver variables over the domain in pc^2, per unspecified out-of-plane depth; the `mass` column is an integral of n, **not grams**. The energy ledger includes integrated net radiative loss; periodic conduction and hydro fluxes conserve their global contribution. Momentum columns use the same solver convention. states.npz contains initial/final conserved states and cell-center coordinates; model.primitive converts these to density, temperature, velocity and pressure.

The Field proxy sqrt[kappa T/(n^2 Lambda)] is an energy-balance estimate, not an exact marginal wavelength or measured interface width. Its minimum is reported at sampled history times, not every stage. The coarse run tests aggregate pressure and cold mass fraction, not pointwise morphology or continuum convergence. Four cells in a proxy is not a proof of convergence. This finite box and enhanced transport intentionally limit interpretation.

The homepage respects reduced-motion preferences, offers play/pause and uses a real-frame poster. Density colors are fixed over the full sequence; spatial enlargement is only display interpolation. Faint lines mark the actual 256 x 160 cell boundaries; slightly stronger guides occur every eight cells. These are display overlays, registered to the computational domain and omitted from the annotation strip. Looping visibly returns to initial time; it is not a periodic physical orbit.

## Attribution

Koyama & Inutsuka (2002), ApJ 564, L97; corrected analytic heating/cooling fit documented in the archive's KI2002 chapter. This particular 2D initial value problem and enlarged conductivity are new illustrative choices for The Art of Thermal Instability, November 2026. Original implementation code is MIT under the project LICENSE. Generated figures and data follow the separate terms in LICENSING.md.
