"""Steady plane-parallel fronts from Iwasaki & Inutsuka (2012), Eqs 5-16.

The classifier follows the outward CNM saddle on the plus-density branch of
Eq. 11. Inputs outside that branch are masked, never assigned by extrapolation.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq

KB=1.380649e-16
MH=1.6735575e-24
R=KB/MH
GAMMA=5/3
CP=GAMMA*R/(GAMMA-1)
HEATING=2e-26
FLUX_UNIT=MH*1e5  # mH * (1 cm^-3) * (1 km/s), g cm^-2 s^-1
F_SCALE=1e-9  # energy-flux scale, erg cm^-2 s^-1; only rescales the ODE


def cooling_ratio(t):
    t=np.asarray(t,dtype=float)
    return 1e7*np.exp(-118400/(t+1000))+.014*np.sqrt(t)*np.exp(-92/t)


def cooling_ratio_derivative(t):
    t=np.asarray(t,dtype=float)
    a=1e7*np.exp(-118400/(t+1000));b=.014*np.sqrt(t)*np.exp(-92/t)
    return a*118400/(t+1000)**2+b*(.5/t+92/t**2)


def kappa(t):return 2500*np.sqrt(t)

def loss(t,n):return HEATING*n*(n*cooling_ratio(t)-1)

def equilibrium_pressure(t):return t/cooling_ratio(t)


def turning_temperatures():
    derivative=lambda t:1/t-float(cooling_ratio_derivative(t)/cooling_ratio(t))
    return brentq(derivative,50,1000),brentq(derivative,1000,10000)


T_COLD_TURN,T_WARM_TURN=turning_temperatures()
P_MIN=float(equilibrium_pressure(T_COLD_TURN))
P_MAX=float(equilibrium_pressure(T_WARM_TURN))


def equilibrium_temperatures(p):
    """Return cold, unstable, warm roots where they exist; pressure is P/kB."""
    fun=lambda t:p*float(cooling_ratio(t))/t-1
    intervals=[(2.,T_COLD_TURN),(T_COLD_TURN,T_WARM_TURN),(T_WARM_TURN,20000.)]
    return [brentq(fun,a,b,xtol=1e-11) for a,b in intervals if fun(a)*fun(b)<0]


def saturation_integral(p):
    tc,_,tw=equilibrium_temperatures(p)
    # Scaling avoids a dimensionful absolute quadrature tolerance at a zero.
    return quad(lambda t:kappa(t)*loss(t,p/t)/F_SCALE**2,tc,tw,epsabs=1e-9,epsrel=1e-9,limit=150)[0]


def saturation_pressure():return brentq(saturation_integral,2000,4000,xtol=1e-7)


@dataclass
class Front:
    p: float
    flux_number: float
    # flux_number = -j / [mH * (1 cm^-3) * (1 km/s)], positive for condensation

    def __post_init__(self):
        if self.p<=P_MIN or self.flux_number<0:raise ValueError('Require CNM pressure > Pmin and nonnegative condensation flux')
        self.tc=equilibrium_temperatures(self.p)[0]
        self.nc=self.p/self.tc;self.rhoc=MH*self.nc;self.j=-self.flux_number*FLUX_UNIT
        self.momentum=self.p*KB+self.j**2/self.rhoc
        self.tcrit=self.momentum**2/(4*R*self.j**2) if self.j else np.inf
        self.initial_isothermal_mach=self.flux_number/(self.nc*np.sqrt(R*self.tc)/1e5)

    def discriminant(self,t):return self.momentum**2-4*R*self.j**2*t

    def density_number(self,t):
        return (self.momentum+np.sqrt(np.maximum(self.discriminant(t),0)))/(2*KB*t)

    def pressure(self,t):return KB*t*self.density_number(t)

    def cooling(self,t):return loss(t,self.density_number(t))

    def cooling_derivative(self,t):
        n=self.density_number(t);dp=-R*self.j**2/np.sqrt(self.discriminant(t))
        dn=dp/(KB*t)-n/t
        return HEATING*((2*n*cooling_ratio(t)-1)*dn+n*n*cooling_ratio_derivative(t))

    def advective_coefficient(self,t):
        # kappa times the coefficient of theta' in equation (14).
        return R*self.j/2*((GAMMA+1)/(GAMMA-1)+self.momentum/np.sqrt(self.discriminant(t)))

    def eigenvalues(self,t):
        a=self.advective_coefficient(t)/kappa(t);b=self.cooling_derivative(t)/kappa(t)
        return np.roots([1,-a,-b])

    def outward_slope(self):
        a=self.advective_coefficient(self.tc)/kappa(self.tc);b=self.cooling_derivative(self.tc)/kappa(self.tc)
        if b<=0 or self.initial_isothermal_mach>=1:return None
        # Stable evaluation of the small positive CNM eigenvalue.
        return 2*b/(-a+np.sqrt(a*a+4*b))

    def thermal_roots(self,maximum_temperature=20000):
        lower=self.tc*(1+1e-6);upper=min(self.tcrit*(1-1e-9),maximum_temperature)
        if upper<=lower:return []
        ts=np.geomspace(lower,upper,1201);values=self.cooling(ts)
        return [brentq(self.cooling,a,b,xtol=1e-10) for a,b,x,y in zip(ts[:-1],ts[1:],values[:-1],values[1:]) if x*y<0]


def phase_flux_squared(front,temperature,forward,midpoint,rtol=3e-8):
    """Integrate F², F=kappa*T', from a saddle into a matching temperature."""
    a=front.advective_coefficient(temperature)/kappa(temperature)
    b=front.cooling_derivative(temperature)/kappa(temperature)
    lam=(a+(1 if forward else -1)*np.sqrt(a*a+4*b))/2
    dt=(1 if forward else -1)*temperature*1e-5
    start=temperature+dt;w0=(lam*kappa(temperature)*dt/F_SCALE)**2
    def rhs(t,w):
        return [(2*front.advective_coefficient(t)*F_SCALE*np.sqrt(max(w[0],0))+2*kappa(t)*front.cooling(t))/F_SCALE**2]
    def stop(t,w):return w[0]
    stop.terminal=True;stop.direction=-1
    sol=solve_ivp(rhs,(start,midpoint),[w0],method='DOP853',events=stop,rtol=rtol,atol=1e-15)
    if not sol.success:raise RuntimeError(sol.message)
    return max(0.,float(sol.y[0,-1]))


def zp_mismatch(p,flux_number,rtol=3e-8):
    front=Front(p,flux_number);roots=front.thermal_roots()
    if len(roots)<2:raise ValueError('No distinct unstable/warm pair for matching')
    tu,tw=roots[:2]
    return phase_flux_squared(front,front.tc,True,tu,rtol)-phase_flux_squared(front,tw,False,tu,rtol)


def zp_flux(p,psat=None,rtol=3e-8):
    psat=saturation_pressure() if psat is None else psat
    if p<=psat:return 0.  # evaporation lies outside this condensation-only plot
    if p>=P_MAX:raise ValueError('ZP pair requires pressure below Pmax')
    hi=.01
    while zp_mismatch(p,hi,rtol)>0:
        hi*=2
        if hi>1:raise RuntimeError('Could not bracket ZP flux')
    return brentq(lambda j:zp_mismatch(p,j,rtol),0,hi,xtol=1e-10)


LABELS={0:'C-U-o',1:'C-U',2:'C-c',3:'C-C',4:'below condensation ZP line',5:'outside outgoing CNM branch',6:'unresolved / outside T window'}


def classify(front,zp_boundary=0,rtol=2e-6,initial_offset=1e-5,temperature_limit=20000):
    """Classify by the global trajectory, with a converged near-root test.

    An unstable-phase node alone does NOT establish C-U (paper footnote 1).
    For real eigenvalues follow log(F) to the first other equilibrium. A
    trajectory captured there has F proportional to its temperature distance;
    an overshoot has finite F. Ambiguous convergence is explicitly unresolved.
    """
    if front.flux_number<=zp_boundary*(1+1e-4) and zp_boundary>0:return 4,dict(reason='Below or on the condensation ZP separatrix')
    lam=front.outward_slope()
    if lam is None:return 5,dict(reason='Eq. 11 plus-density branch lacks an outgoing CNM saddle',initial_isothermal_mach=float(front.initial_isothermal_mach))
    roots=front.thermal_roots(maximum_temperature=temperature_limit)
    if not roots:
        if front.tcrit<temperature_limit:return 2,dict(reason='No intervening thermal fixed point before the critical temperature',T_terminal_K=float(front.tcrit))
        return 6,dict(reason=f'No terminal point within the audited {temperature_limit:g} K integration window')
    target=roots[0];eigen=front.eigenvalues(target)
    if np.any(abs(np.imag(eigen))>1e-30) and target>T_COLD_TURN:
        return 0,dict(reason='Unstable-phase spiral; temperature overshoot',T_terminal_K=float(target))
    dt=min(front.tc*initial_offset,(target-front.tc)*1e-3);start=front.tc+dt
    Finit=lam*kappa(front.tc)*dt
    def rhs(t,y):
        F=F_SCALE*np.exp(np.clip(y[0],-200,200))
        return [front.advective_coefficient(t)/F+kappa(t)*front.cooling(t)/F**2]
    def jac(t,y):
        F=F_SCALE*np.exp(np.clip(y[0],-200,200))
        return np.array([[-front.advective_coefficient(t)/F-2*kappa(t)*front.cooling(t)/F**2]])
    last_t=start;last_y=np.log(Finit/F_SCALE);values=[];previous_kind=None
    for eps in [1e-3,1e-4,1e-5,1e-6,1e-7]:
        end=target-eps*(target-front.tc)
        sol=solve_ivp(rhs,(last_t,end),[last_y],method='Radau',jac=jac,rtol=rtol,atol=rtol,max_step=(target-front.tc)/8)
        if not sol.success:return 6,dict(reason='ODE failed: '+sol.message)
        last_t=end;last_y=sol.y[0,-1];F=F_SCALE*np.exp(last_y);values.append(F)
        if len(values)>=2:
            exponent=np.log(values[-1]/values[-2])/np.log(.1)
            kind='capture' if .8<exponent<1.2 else 'cross' if abs(exponent)<.1 else None
            confirmed=kind is not None and kind==previous_kind
            previous_kind=kind
            if confirmed and kind=='capture':
                return (3 if target<T_COLD_TURN else 1),dict(reason='Converged approach to thermal fixed point',T_terminal_K=float(target),flux_distance_exponent=float(exponent),distance_fraction=eps)
            if confirmed and kind=='cross':
                if target<T_COLD_TURN:return 6,dict(reason='Overshoot of a cold fixed point requires further trajectory tracing')
                return 0,dict(reason='Finite heat flux across unstable fixed-point temperature',T_terminal_K=float(target),flux_distance_exponent=float(exponent),distance_fraction=eps)
    return 6,dict(reason='Near-root asymptotics did not separate capture from overshoot',T_terminal_K=float(target))
