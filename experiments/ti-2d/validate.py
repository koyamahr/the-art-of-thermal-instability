"""Independent continuum benchmarks; run with Python 3.9+, NumPy and SciPy."""
import json,tempfile,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from model import *

def validate(binary,c,root):
    checks=[]
    def check(name,value,limit):
        checks.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(value<limit)))
    cfg=dict(c,domain_pc=[4.8,.3]);nx=96;ny=4;x=(np.arange(nx)+.5)*4.8/nx;k=2*np.pi/4.8
    n=np.ones((ny,nx));wave=np.cos(k*x)[None,:];t0=1000.;p0=B*t0
    def go(name,u,mode,end,**kw):return run(binary,u,cfg,root/name,end=end,mode=mode,frames=8,**kw)
    def amp(a):return 2*np.mean(a*np.exp(-1j*k*x)[None,:])
    # Advected entropy/contact wave, analytic translation.
    u=conservative(n*(1+.1*wave),t0/(1+.1*wave),.5,0)
    f,h,_=go('contact',u,1,1.)
    check('advected contact Fourier error',abs(amp(f[...,0])/(.1*np.exp(-1j*k*.5))-1),.025)
    check('contact pressure preservation',np.max(abs(primitive(f)[4]/p0-1)),1e-10)
    # Constant-coefficient heat equation: continuum exponential decay.
    kap=3e6;chi=KUNIT*kap*(GAMMA-1)/B;end=.2/(chi*k*k)
    f,h,_=go('diffusion',conservative(n,t0+wave),2,end,base=kap,cold=0)
    check('heat equation Fourier error',abs(amp(primitive(f)[1])/(np.exp(-.2))-1),.001)
    check('conduction energy conservation',abs(h[-1]['total_energy']/h[0]['total_energy']-1),1e-11)
    # Independent adaptive ODE, no spatial transport.
    nn=2.8*n;end=.2
    exact=solve_ivp(lambda _,t:-(GAMMA-1)*HEAT/B*(2.8*ratio(t)-1),[0,end],[2000.],rtol=1e-11,atol=1e-10).y[0,-1]
    f,h,_=go('cooling',conservative(nn,2000.),3,end)
    check('isochoric cooling versus adaptive ODE',abs(np.mean(primitive(f)[1])/exact-1),.0005)
    # Linear compressible TI eigenmode: independent Fourier-space 3x3 matrix.
    n0=2.8;t=brentq(lambda t:n0*ratio(t)-1,300,4000);p=B*n0*t
    d=HEAT*n0*n0*derivative(t)+KUNIT*kappa(t,c)*k*k
    qn=HEAT*(2*n0*ratio(t)-1)
    a=np.array([[0,-1j*k*n0,0],[0,0,-1j*k/n0],[(GAMMA-1)*(d*t/n0-qn),-1j*k*GAMMA*p,-(GAMMA-1)*d/(B*n0)]],complex)
    vals,vec=np.linalg.eig(a);idx=np.argmax(vals.real);s=vals[idx];v=vec[:,idx];v*=n0*1e-5/v[0]
    z=np.exp(1j*k*x)[None,:];dn=(v[0]*z).real;vel=(v[1]*z).real;dp=(v[2]*z).real
    u=conservative(n*n0+dn,(p+dp)/(B*(n*n0+dn)),vel,0);end=1.
    f,h,_=go('linear-ti',u,0,end);nf,tf,vf,_,pf=primitive(f)
    for name,field,expected in [('density',nf,v[0]),('velocity',vf,v[1]),('pressure',pf,v[2])]:
        check('linear TI '+name+' eigenmode error',abs(amp(field)/(expected*np.exp(s*end))-1),.025)
    check('positive TI eigenvalue',0 if s.real>0 else 1,.5)
    check('radiative energy ledger',abs((h[-1]['total_energy']-h[0]['total_energy']+h[-1]['integrated_net_radiative_loss'])/h[0]['total_energy']),1e-10)
    return checks

if __name__=='__main__':
    c=json.loads((HERE/'config.json').read_text())
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp);checks=validate(compile_solver(root),c,root)
    report=dict(checks=checks,code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'solver.cpp',HERE/'model.py',HERE/'validate.py']},config_sha256=hashlib.sha256((HERE/'config.json').read_bytes()).hexdigest())
    if '--check-only' not in __import__('sys').argv: (HERE/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(checks,indent=2));assert all(c['passed'] for c in checks)
