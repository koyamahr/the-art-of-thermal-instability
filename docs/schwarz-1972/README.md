# Schwarz, McCray & Stein (1972): an early-time linear prototype

Target: *Formation of Clouds in a Cooling Interstellar Medium*, ApJ175,
673–686. Figure2 is the first calculation target. Nonlinear Figures3–6
remain a separate, unfinished stage.

## Run

Keep reproduce.py, config.json, cooling-readings.csv, reference-readings.csv,
and requirements.txt together. The source-audit.json documents the scope.

```
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduce.py
```

Outputs: linear.png/svg/csv, cooling.csv, reference-comparison.csv,
and validation.json under output/. No PDF or network access is needed
once dependencies are installed. Public source rasters are not redistributed.

## Equations and assumptions

The background obeys Eqs16–18, and four perturbations per wavelength obey
Eqs19–24: density, pressure, scaled velocity divergence, and ionization.
Temperature perturbations are reconstructed from the equation of state.
The code includes neutral thermal conduction K=2500sqrt(T) in cgs units,
H recombination and the changing particle count in the energy equation.
There is no continuing heating or cosmic-ray ionization.

Initial nH=.3cm^-3, T=8200K, x=.05, neutral He/H=.1, gamma=5/3;
wavelengths=.3,1,3,10pc. Density amplitudes are normalized to1. Pressure
perturbations, velocity, and ionization perturbations start at zero.
This selects an initially isobaric, uniformly ionized perturbation. The
nonlinear section explicitly states an isobaric initial perturbation, but
its exact use for Figure2 has not been established. No time rescaling or
coefficient fitting is used to match the target curves.

The ionization coefficient uses the excited-state hydrogen approximation
from Hjellming1966 Eqs56–57. The precise numerical recombination fit used
in the target paper is not specified there and remains unverified.

## Cooling input and stopping rule

Cooling is reconstructed from39 readings of the three Figure1 curves,
x=.001,.01,.1. It is **not** an independent atomic-rate reconstruction.
PCHIP interpolates logLambda(logT); Lambda is piecewise linear in x,
motivated by the electron/neutral collider contributions for fixed trace
ions. The logarithmic x interpolation alternative changes the stopping
time from about1.679 to1.509Myr. This is a sensitivity probe, not a
statistical uncertainty interval or a fitted correction.

The accepted output stops at T=13.929K, the lowest temperature with an
input reading, after about1.679Myr. No colder input is extrapolated.
This excludes the late-time portion of Figure2, including the later
behavior of the3pc and10pc modes. Maxima reported at the final sample
are maxima **within the computed interval**, not established mode peaks.
Negative linear density amplitudes remain in the CSV and are omitted
on the log plot; they do not mean that the physical density is negative.

Raw pixel coordinates, source hash, crop origins, calibration ticks,
interpolation and practical reading envelopes are retained. Twenty-four
sparse Figure2 readings provide a comparison; points after the domain
stop have no model prediction and are unscored. No independent second
reading has been performed. RMS discrepancies are descriptive and are
not numerical pass tolerances.

## Checks and historical caveats

The31 checks include axis calibration, interpolation of the source points,
analytic coefficient derivatives, Eq19 versus a centered-difference
Jacobian of the original primitive pressure equation, independent DOP853
and Radau integrations, background domain/monotonicity, and conduction
controls. No assertion of full historical-figure agreement follows.

The SMS Fe abundance is5e-5 perH, unlike the8e-6 used in the current FGH
prototype. FGH and KI rates are not silently substituted. Jura & Dalgarno
1972 Table1 supplies several electron-impact terms but does not alone
supply all neutral O/H channels needed here. Recovering every rate is a
separate task.

Eq7 uses radiated power per gram nLambda/(1.4mH), while the Figure1
caption omits1.4. This run follows Eqs7–9 without correcting the plotted
normalization to fit Figure2. The ambiguity is recorded. The note added
in proof (p686) revises the interpretation of trace-element ionization
and cooling; it does not supply a replacement numerical rate table.

## Next stage

Recover or extend the cold input with explicit source justification,
confirm the recombination coefficient and initial linear mode, then
revisit the full Figure2. For Figures3–6, implement periodic1D flow with
a3pc wavelength and10% initially isobaric density modulation, including
conduction, recombination, shock handling and grid/time convergence.
A large linear amplification alone cannot reproduce the finite-amplitude
cloud contrast, thickness or shock.
