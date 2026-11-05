# Koyama & Inutsuka (2002): corrected equilibrium benchmark

This is a new equilibrium diagnostic of the corrected analytic fit, not a rerun
of the 2002 two-dimensional simulation or a digitization of an original figure.

Run with Python 3.9 or newer:

```sh
python -m pip install -r requirements.txt
python reproduce.py --output output
```

Download equilibrium.py, reproduce.py, config.json, source-audit.json and
requirements.txt into one directory. No source PDFs are required to run.
The repository's publish.py also cross-checks the actual Iwasaki front code;
the standalone run checks frozen values with the reference code hash recorded.

## Model and correction

Use n=rho/m_H (pure atomic H), P/k_B=nT and net volumetric cooling
Q=n^2 Lambda(T)-n Gamma, positive for cooling. Gamma=2e-26 erg per H per second.
The corrected ratio Lambda/Gamma in cm^3 is

    1e7 exp[-118400/(T+1000)] + 0.014 sqrt(T) exp[-92/T].

The inspected arXiv v1 Eq.4 instead has 114800 and 14. Vazquez-Semadeni et al.
(2007), ApJ657, p.873, Eq.4 and footnote5, identify two typos and attribute the
correction to H.Koyama. This is a correction-provenance citation, not an addition
of that simulation to the atlas. No separate formal journal erratum was located;
this does not establish its absence. Publisher PDF access failed; we do not claim
to have independently inspected its typesetting. See source-audit.json for URLs,
editions, hashes and the visually checked Iwasaki2012 Eqs4--7.

## Stability and ZP selection

Equilibrium: n=Gamma/Lambda, P/k_B=T Gamma/Lambda. Set beta=d ln Lambda/d ln T.
At fixed P, (dQ/dT)_P=n Gamma (beta-1)/T at equilibrium. Thus the local
isobaric mode is stable for beta>1, unstable for beta<1 and marginal for beta=1.
The unstable interval is 184--5039.419 K. This label omits conduction's
wavelength dependence. Isochoric loss derivative is positive throughout the
plotted interval; isobaric and isochoric criteria must not be interchanged.
The displayed 10--10000 K range is a plotting choice, not a validation of the
fit's chemical accuracy over that interval.

Uniform equilibria have zero temperature gradient, so conduction does not move
the curve. For a planar static front, d(kappa dT/dx)/dx=Q and zero conductive
flux at both endpoints instead require int[Tc,Tw] kappa(T) Q(T,P) dT=0.
This is not an unweighted equal-area construction on the P--n plot.
We use kappa=2500 sqrt(T) erg cm^-1 K^-1 s^-1, as in Iwasaki2012.

P_sat/k_B=2822.884297 K cm^-3:

- Cold endpoint: T=54.980188 K, n=51.343664 cm^-3.
- Unstable intersection: T=1364.881273 K, n=2.068227 cm^-3.
- Warm endpoint: T=6417.132178 K, n=0.439898 cm^-3.

The intermediate intersection is not a third asymptotic endpoint of the same
cold-to-warm static front. The common-pressure equilibrium range is
1597.453101 < P/k_B < 5007.000118; a static planar conductive front selects
one pressure within it. Iwasaki prints 1596 (text), 1598 (Figure1), 5012 and
2823. Small residuals are retained in validation.json; no coefficients fitted.

Outputs: equilibrium.csv, saturation.csv, coefficient-comparison.csv, SVG/PNG
figures and validation.json. Checks cover energy closure, analytic vs finite
difference derivatives, independent equilibrium inversion, parcel response,
independent conductive quadrature and agreement with the Iwasaki code.


## Pressure-dependent static-front potential

Run `python potential.py --output output` with equilibrium.py, config-potential.json and the same pinned dependencies. No multidimensional calculation is required. potential.svg/png compares 0.9, 1 and 1.1 Psat in log temperature and in the natural conduction coordinate; potential-mobile.svg/png stacks temperature panels for narrow displays. potential.csv preserves every curve, while validation-potential.json records all equilibria, densities, well depths, normalization, independent checks and source/configuration hashes.

The coordinate is theta=(5000/3)T^(3/2), in erg cm^-1 s^-1. V(theta(T),P)=integral_Tcold^T kappa Q dT, in erg^2 cm^-4 s^-2. Vcold=0 independently at each pressure. All panels divide by the same saturation barrier, 3.425672655e-18 erg^2 cm^-4 s^-2. This effective potential is not thermodynamic free energy. A static zero-flux connection requires Vwarm=Vcold; dV/dtheta=Q and equilibrium curvature is Q_T/kappa. A logarithmic-temperature axis changes well widths but preserves extrema and depth differences.

Independent Simpson integration and 16001-to-32001-point refinement check the adaptive quadrature; finite differences check the derivative in theta and all equilibrium curvatures. Below saturation the warm minimum is deeper; above it the cold minimum is deeper. The report retains the small signed cancellation residual at Psat. No time-dependent compressible front velocity is computed or inferred.
