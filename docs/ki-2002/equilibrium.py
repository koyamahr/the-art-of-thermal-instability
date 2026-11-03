"""Corrected KI (2002) fit; pure atomic H and constant heating per H atom.

Temperatures in kelvin; density n=rho/m_H in cm^-3; pressure argument is P/kB
in K cm^-3. Net loss is positive for cooling, in erg cm^-3 s^-1.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

HEATING = 2e-26
KB = 1.380649e-16


def ratio(t, printed=False):
    """Lambda/Gamma in cm^3. printed=True is the uncorrected arXiv v1 Eq.4."""
    t = np.asarray(t, dtype=float)
    if np.any(t <= 0):
        raise ValueError('Temperature must be positive')
    a, b = (114800., 14.) if printed else (118400., .014)
    return 1e7*np.exp(-a/(t+1000)) + b*np.sqrt(t)*np.exp(-92/t)


def slope(t):
    """d ln Lambda / d ln T for the corrected fit."""
    t = np.asarray(t, dtype=float)
    a = 1e7*np.exp(-118400/(t+1000))
    b = .014*np.sqrt(t)*np.exp(-92/t)
    return (a*118400*t/(t+1000)**2 + b*(.5+92/t))/(a+b)


def loss(t, n):
    return HEATING*n*(n*ratio(t)-1)


def pressure(t):
    return t/ratio(t)


def turns():
    return [brentq(lambda t: slope(t)-1, 50, 1000),
            brentq(lambda t: slope(t)-1, 1000, 10000)]


def roots(p):
    tc, tw = turns()
    intervals = [(10., tc), (tc, tw), (tw, 10000.)]
    f = lambda t: p*ratio(t)/t-1
    return [brentq(f, a, b, xtol=1e-10) for a, b in intervals if f(a)*f(b) < 0]


def static_integral(p):
    """Scaled integral int kappa(T) Q(T,P) dT; zero selects a static front."""
    ts = roots(p)
    if len(ts) != 3:
        raise ValueError('Pressure must have three equilibria inside 10--10000 K')
    return quad(lambda t: 2500*np.sqrt(t)*loss(t, p/t)/1e-18,
                ts[0], ts[-1], epsabs=1e-10, epsrel=1e-10, limit=150)[0]


def saturation_pressure():
    return brentq(static_integral, 2000, 4000, xtol=1e-8)
