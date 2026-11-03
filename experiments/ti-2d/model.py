"""Units, reproducible seeds, and independent continuum diagnostics for TI2D."""
from pathlib import Path
import json, subprocess, time
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

KB=1.380649e-16; MH=1.6735575e-24; PC=3.0856775814913673e18; MYR=31557600e6
MASS=1.37*MH; GAMMA=5/3; VUNIT=PC/MYR; B=KB/(MASS*VUNIT**2)
HEAT=2e-26*MYR/(MASS*VUNIT**2)
KUNIT=MYR/(PC**2*MASS*VUNIT**2)
HERE=Path(__file__).resolve().parent

def ratio(t):
    t=np.asarray(t);return 1e7*np.exp(-118400/(t+1000))+.014*np.sqrt(t)*np.exp(-92/t)

def derivative(t):
    a=1e7*np.exp(-118400/(t+1000));b=.014*np.sqrt(t)*np.exp(-92/t)
    return a*118400/(t+1000)**2+b*(.5/t+92/t**2)

def kappa(t,c):return c['conductivity']['base_cgs']+c['conductivity']['cold_cgs']*(c['conductivity']['reference_K']/t)**2

def primitive(u):
    n=u[...,0];vx=u[...,1]/n;vy=u[...,2]/n
    p=(GAMMA-1)*(u[...,3]-.5*n*(vx*vx+vy*vy));return n,p/(B*n),vx,vy,p

def conservative(n,t,vx=0,vy=0):
    n,t,vx,vy=np.broadcast_arrays(n,t,vx,vy);p=B*n*t
    return np.stack([n,n*vx,n*vy,p/(GAMMA-1)+.5*n*(vx*vx+vy*vy)],axis=-1)

def initial(c,nx=None,ny=None):
    nx=nx or c['grid'][0];ny=ny or c['grid'][1];lx,ly=c['domain_pc'];xx,yy=np.meshgrid((np.arange(nx)+.5)*lx/nx,(np.arange(ny)+.5)*ly/ny)
    rng=np.random.default_rng(c['seed']);pattern=np.zeros_like(xx)
    for i in range(8):
        for j in range(-5,6):
            if i==0 and j<=0:continue
            q=np.sqrt((i/lx)**2+(j/ly)**2)
            if q>1.5:continue
            pattern+=rng.normal()/q**1.5*np.cos(2*np.pi*(i*xx/lx+j*yy/ly)+rng.uniform(0,2*np.pi))
    pattern-=pattern.mean();pattern/=np.std(pattern)
    n=c['mean_n_cm3']*(1+c['density_rms_fraction']*pattern)
    t0=brentq(lambda t:c['mean_n_cm3']*ratio(t)-1,300,4000)
    t=t0*c['mean_n_cm3']/n
    return conservative(n,t)

def compile_solver(directory):
    binary=Path(directory)/'ti2d-solver';binary.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(['c++','-O3','-std=c++17',str(HERE/'solver.cpp'),'-o',str(binary)],check=True)
    return binary

def run(binary,u,c,out,end=None,mode=0,frames=2,base=None,cold=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);u.astype('<f8').tofile(out/'initial.bin')
    ny,nx=u.shape[:2];lx,ly=c['domain_pc'];start=time.perf_counter()
    args=[str(binary),str(nx),str(ny),str(lx),str(ly),str(end if end is not None else c['end_time_Myr']),str(mode),str(c['conductivity']['base_cgs'] if base is None else base),str(c['conductivity']['cold_cgs'] if cold is None else cold),str(c['conductivity']['reference_K']),str(out/'initial.bin'),str(out),str(frames)]
    subprocess.run(args,check=True)
    elapsed=time.perf_counter()-start
    final=np.fromfile(out/'final.bin',dtype='<f8').reshape(ny,nx,4)
    hist=np.genfromtxt(out/'history.csv',delimiter=',',names=True)
    return final,hist,elapsed

def landmarks(c):
    tc=brentq(lambda t:t*derivative(t)/ratio(t)-1,50,1000);tw=brentq(lambda t:t*derivative(t)/ratio(t)-1,1000,10000)
    pm=tc/ratio(tc);px=tw/ratio(tw)
    def roots(p):return [brentq(lambda t:p*ratio(t)/t-1,a,b) for a,b in [(10,tc),(tc,tw),(tw,10000)]]
    def integral(p):
        a,_,b=roots(p)
        return quad(lambda t:kappa(t,c)/c['conductivity']['base_cgs']*(p/t)*(p*ratio(t)/t-1),a,b,epsabs=1e-6,epsrel=1e-9,limit=150)[0]
    ps=brentq(integral,pm*(1+1e-6),px*(1-1e-6))
    return dict(P_sat_over_kB_K_cm3=ps,states=[dict(T_K=t,n_cm3=ps/t) for t in roots(ps)],conductivity=c['conductivity'])
