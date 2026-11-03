"""Schwarz, McCray & Stein1972: early-time Figure2 linear prototype.

Input cooling comes from39 graph readings, not independent atomic-rate recovery.
The run stops at the input temperature/ionization boundary. No extrapolation.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
KB=1.380649e-16
MH=1.6735575e-24
PC=3.0856775814913673e18
MYR=365.25*86400*1e6


def axis_fit(ticks):
    p,v=np.asarray(ticks,float).T;slope,offset=np.polyfit(p,v,1)
    return slope,offset,float(np.max(abs(slope*p+offset-v)))


class Cooling:
    def __init__(self,config):
        d=config['cooling_digitization'];sx,bx,_=axis_fit(d['x_ticks']);sy,by,_=axis_fit(d['y_ticks'])
        rows=list(csv.DictReader((ROOT/'cooling-readings.csv').open()));self.rows=[];self.curves=[]
        self.mode=config.get('x_interpolation','linear_lambda');self.xs=sorted(set(float(r['electron_fraction']) for r in rows))
        for x in self.xs:
            r=[z for z in rows if float(z['electron_fraction'])==x]
            t=10**np.array([sx*float(z['x_px'])+bx for z in r]);lam=10**np.array([sy*float(z['y_px'])+by for z in r])
            self.curves.append(PchipInterpolator(np.log(t),np.log(lam),extrapolate=False))
            self.rows.extend(zip(np.full(len(t),x),t,lam))
        self.tmin=float(min(t));self.tmax=float(max(t));self.logxs=np.log(self.xs)
    def evaluate(self,t,x):
        # A solver can trial outside an event surface. Clamp coefficient evaluations;
        # the accepted solution is truncated at the source boundary event.
        lt=np.log(np.clip(t,self.tmin,self.tmax));lx=np.log(np.clip(x,self.xs[0],self.xs[-1]))
        i=int(np.clip(np.searchsorted(self.logxs,lx)-1,0,len(self.xs)-2));a,b=self.logxs[i:i+2]
        w=(lx-a)/(b-a);y0=float(self.curves[i](lt));y1=float(self.curves[i+1](lt))
        st0=float(self.curves[i].derivative()(lt));st1=float(self.curves[i+1].derivative()(lt))
        if self.mode=='log_lambda':
            return float(np.exp((1-w)*y0+w*y1)),(1-w)*st0+w*st1,(y1-y0)/(b-a)
        xclip=np.exp(lx);xa,xb=self.xs[i:i+2];w=(xclip-xa)/(xb-xa);l0,l1=np.exp([y0,y1]);value=(1-w)*l0+w*l1
        return float(value),float(((1-w)*l0*st0+w*l1*st1)/value),float(xclip*(l1-l0)/((xb-xa)*value))


def recombination(t):
    lam=157890/np.asarray(t,float)
    f=.4288+.5*np.log(lam)+.469*lam**(-1/3)-(.7776-.4896/lam+.4796/lam**2-.057/lam**3-.405/lam**4)
    df=.5/lam-.469/3*lam**(-4/3)-(.4896/lam**2-.9592/lam**3+.171/lam**4+1.62/lam**5)
    return 5.197e-14*np.sqrt(lam)*f,-.5-lam*df/f


def rhs(time,y,config,cooling,conduction=True):
    initial=config['initial'];n0=initial['nH_cm3'];h=1+initial['helium_per_H'];g=initial['gamma']
    t,x=np.exp(y[:2]);lam,st,sx=cooling.evaluate(t,x);alpha,sa=recombination(t)
    rc=MYR*n0*(g-1)*lam/((h+x)*KB*t);rr=MYR*n0*x*alpha;r=x/(h+x)
    result=np.zeros_like(y);result[0]=-rc+(2-g)*r*rr;result[1]=-rr
    for i,wave in enumerate(config['wavelength_pc']):
        n,p,q,z=y[2+4*i:6+4*i];temp=p-n-r*z;k=2*np.pi/(wave*PC)
        sound2=g*(h+x)*KB*t/(1.4*MH)
        dc=MYR*config['conductivity_cgs_prefactor']*np.sqrt(t)*k*k*(g-1)/(n0*(h+x)*KB) if conduction else 0.
        result[2+4*i]=-q
        result[3+4*i]=-rc*((2-st)*(n-p)+p+(sx-r*st)*z)-(g-1)*r*rr*((2-r*(1+sa))*z+(1-sa)*n+sa*p)-g*q-dc*temp
        result[4+4*i]=MYR**2*k*k*sound2/g*p
        result[5+4*i]=-rr*(n+z+sa*temp)
    return result


def integrate(config,cooling,method='DOP853',rtol=2e-9,conduction=True):
    a=config['initial'];y=[np.log(a['T_K']),np.log(a['electron_fraction'])]
    for w in config['wavelength_pc']:y.extend([a['relative_density'],a['relative_pressure'],a['scaled_velocity_divergence'],a['relative_ionization']])
    def low_temperature(t,y):return y[0]-np.log(cooling.tmin)
    def low_ionization(t,y):return y[1]-np.log(cooling.xs[0])
    low_temperature.terminal=True;low_temperature.direction=-1
    low_ionization.terminal=True;low_ionization.direction=-1
    sol=solve_ivp(lambda t,y:rhs(t,y,config,cooling,conduction),[0,config['time_max_Myr']],y,method=method,rtol=rtol,atol=rtol*1e-3,dense_output=True,events=[low_temperature,low_ionization])
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def primitive_pressure_rate(y,config,cooling,amp=0.,mode=None,wave=.3):
    """Original pressure, continuity, ionization equations for Jacobian checking.

    Physical local source variations are centered-differenced, with the moving
    background normalization subtracted. This does not reuse Eq19.
    """
    init=config['initial'];h=1+init['helium_per_H'];g=init['gamma'];n0=init['nH_cm3'];t0,x0=np.exp(y[:2]);r=x0/(h+x0)
    nv,pv,qv,zv=mode if mode is not None else [0,0,0,0]
    temp=pv-nv-r*zv;n=n0*(1+amp*nv);t=t0*(1+amp*temp);x=x0*(1+amp*zv)
    pressure_ratio=(n/n0)*(t/t0)*(h+x)/(h+x0)
    lam,_,_=cooling.evaluate(t,x);alpha,_=recombination(t)
    cool=MYR*(g-1)*n*n*lam/(n0*(h+x0)*KB*t0)
    rec=MYR*(g-1)*n*n*x*x*alpha*t/(n0*(h+x0)*t0)
    rate=-g*pressure_ratio*amp*qv-cool-rec
    k=2*np.pi/(wave*PC);dc=MYR*config['conductivity_cgs_prefactor']*np.sqrt(t0)*k*k*(g-1)/(n0*(h+x0)*KB)
    return rate-dc*amp*temp


def main():
    c=json.loads((ROOT/'config.json').read_text());cool=Cooling(c);checks=[];out=ROOT/'output';out.mkdir(exist_ok=True)
    def check(name,error,tol):checks.append(dict(name=name,error=float(error),tolerance=tol,passed=bool(error<tol)))
    def condition(name,v):checks.append(dict(name=name,passed=bool(v)))
    for axis in ['x','y']:check('Figure1 '+axis+' calibration residual in dex',axis_fit(c['cooling_digitization'][axis+'_ticks'])[2],.03)
    for axis in ['x','y']:check('Figure2 '+axis+' calibration residual (Myr or dex)',axis_fit(c['growth_digitization'][axis+'_ticks'])[2],.03)
    condition('Cooling source contains three13-point curves',len(cool.rows)==39 and len(cool.xs)==3)
    errors=[abs(cool.evaluate(t,x)[0]/v-1) for x,t,v in cool.rows]
    check('PCHIP interpolates every supplied graph reading',max(errors),1e-12)
    for temp,x in [(50,.003),(500,.02),(7800,.05)]:
        value,st,sx=cool.evaluate(temp,x);eps=1e-5
        check('Cooling T derivative at '+str(temp),abs((np.log(cool.evaluate(temp*np.exp(eps),x)[0])-np.log(cool.evaluate(temp*np.exp(-eps),x)[0]))/(2*eps)-st),1e-7)
        check('Cooling x derivative at '+str(temp),abs((np.log(cool.evaluate(temp,x*np.exp(eps))[0])-np.log(cool.evaluate(temp,x*np.exp(-eps))[0]))/(2*eps)-sx),1e-7)
        alpha,sa=recombination(temp)
        check('Recombination slope at '+str(temp),abs((np.log(recombination(temp*np.exp(eps))[0])-np.log(recombination(temp*np.exp(-eps))[0]))/(2*eps)-sa),1e-8)
        y=np.array([np.log(temp),np.log(x),.4,.2,-.03,.1]);cc=dict(c,wavelength_pc=[.3])
        eps=1e-5;mode=y[2:];zero=primitive_pressure_rate(y,cc,cool)
        diff=(primitive_pressure_rate(y,cc,cool,eps,mode)-primitive_pressure_rate(y,cc,cool,-eps,mode))/(2*eps)-mode[1]*zero
        target=rhs(0,y,cc,cool)[3]
        check('Eq19 agrees with primitive-equation Jacobian at '+str(temp),abs(diff-target)/max(1.,abs(target)),2e-8)
    sol=integrate(c,cool);reference=integrate(c,cool,method='Radau',rtol=5e-10)
    check('Independent integrators agree on stopping time',abs(sol.t[-1]-reference.t[-1]),2e-6)
    times=np.linspace(0,min(sol.t[-1],reference.t[-1]),c['time_samples']);y=sol.sol(times);ref=reference.sol(times)
    check('Background temperature agrees between integrators',np.max(abs(np.exp(y[0]-ref[0])-1)),2e-6)
    check('Background ionization agrees between integrators',np.max(abs(np.exp(y[1]-ref[1])-1)),2e-6)
    condition('Accepted output remains inside cooling-temperature domain',np.min(np.exp(y[0]))>=cool.tmin*(1-1e-8))
    condition('Accepted output remains inside ionization domain',np.min(np.exp(y[1]))>=cool.xs[0]*(1-1e-8))
    condition('Background cools and recombines monotonically',np.all(np.diff(y[0])<0) and np.all(np.diff(y[1])<0))
    no_conduction=integrate(c,cool,conduction=False);yn=no_conduction.sol(times)
    d=c['growth_digitization'];sx,bx,_=axis_fit(d['x_ticks']);sy,by,_=axis_fit(d['y_ticks'])
    readings=[]
    for r in csv.DictReader((ROOT/'reference-readings.csv').open()):
        time=sx*float(r['x_px'])+bx;value=10**(sy*float(r['y_px'])+by+(1 if r['quantity']=='T' else 0))
        wave=float(r['wavelength_pc']) if r['wavelength_pc'] else None
        model=None
        if 0<=time<=sol.t[-1]:
            state=sol.sol(time);model=float(np.exp(state[0])) if wave is None else float(state[2+4*c['wavelength_pc'].index(wave)])
        readings.append(dict(quantity=r['quantity'],wavelength_pc=wave,time_Myr=float(time),reference_value=float(value),model_value=model,residual_dex=float(np.log10(model/value)) if model is not None and model>0 else None))
    alternative_cooling=Cooling(dict(c,x_interpolation='log_lambda'));alternative=integrate(c,alternative_cooling)
    interpolation_probe=dict(alternative='logLambda linear in logx; a sensitivity alternative, not an uncertainty bound',stopping_time_Myr=float(alternative.t[-1]))
    summary=[]
    fig,axes=plt.subplots(1,2,figsize=(12,5.7),layout='constrained')
    ax=axes[0];ax2=ax.twinx()
    for i,w in enumerate(c['wavelength_pc']):
        v=y[2+4*i];vv=ref[2+4*i]
        check(str(w)+'pc density agrees between integrators',np.max(abs(v-vv)/(1+abs(vv))),2e-5)
        line,=ax.semilogy(times,np.where(v>0,v,np.nan),label=str(w)+' pc',lw=1.8)
        pp=[r for r in readings if r['wavelength_pc']==w]
        ax.plot([r['time_Myr'] for r in pp],[r['reference_value'] for r in pp],'o',mfc='white',mec=line.get_color(),ms=4)
        j=int(np.argmax(v));summary.append(dict(wavelength_pc=w,max_amplification_in_window=float(v[j]),peak_time_Myr=float(times[j]),minimum_amplification=float(np.min(v)),end_amplification=float(v[-1])))
    condition('0.3pc density initially decreases with conduction',np.min(y[2,:max(2,len(times)//4)])<.95)
    condition('Removing conduction changes the short-wave solution',np.max(abs(yn[2]-y[2]))>.1)
    check('Conduction does not change the uniform cooling background',np.max(abs(yn[:2]-y[:2])),2e-6)
    ax2.semilogy(times,np.exp(y[0]),'k--',lw=1,label='Background T')
    pp=[r for r in readings if r['quantity']=='T'];ax2.plot([r['time_Myr'] for r in pp],[r['reference_value'] for r in pp],'kx',ms=4)
    ax.axvspan(sol.t[-1],2.1,color='gray',alpha=.12);ax.set_xlim(0,2.1)
    ax.set(xlabel='Time (Myr)',ylabel=r'$n_1(t)/n_1(0)$',title='Figure2: early-time linear prototype',ylim=(.1,300));ax2.set(ylabel='Background temperature (K)',ylim=(10,1e4));ax.legend(frameon=False,loc='lower left',fontsize=8);ax.grid(alpha=.15)
    for x in cool.xs:
        pts=[(t,v) for xx,t,v in cool.rows if xx==x];t=np.geomspace(cool.tmin,cool.tmax,300)
        line,=axes[1].loglog(t,[cool.evaluate(tt,x)[0] for tt in t],label='x='+str(x))
        axes[1].loglog(*np.array(pts).T,'o',color=line.get_color(),ms=3)
    axes[1].set(xlabel='Temperature (K)',ylabel=r'$\Lambda$ (erg cm$^3$ s$^{-1}$)',title='Input:39 readings from published Figure1');axes[1].legend(frameon=False);axes[1].grid(alpha=.15)
    fig.suptitle('Schwarz, McCray & Stein (1972): cooling, recombination, and conduction')
    fig.supxlabel('Graph-reconstructed cooling; declared hydrogen recombination approximation; no fitted time rescaling.\nOpen symbols: sparse Figure2 readings; gray: beyond source-data domain. Late-time Figure2 and nonlinear Figures3–6 are not reproduced.',fontsize=9)
    matplotlib.rcParams['svg.hashsalt']='sms72-linear'
    fig.savefig(out/'linear.png',dpi=180);fig.savefig(out/'linear.svg',metadata={'Date':None});plt.close(fig)
    p=out/'linear.svg';p.write_text('\n'.join(z.rstrip() for z in p.read_text().splitlines())+'\n')
    with (out/'linear.csv').open('w',newline='') as f:
        wr=csv.writer(f,lineterminator='\n');wr.writerow(['time_Myr','T0_K','x0']+[q+'_'+str(w)+'pc' for w in c['wavelength_pc'] for q in ['n1','p1','Myr_div_v1','x1']]);wr.writerows(np.column_stack([times,np.exp(y[0]),np.exp(y[1]),*y[2:]]))
    with (out/'cooling.csv').open('w',newline='') as f:
        wr=csv.writer(f,lineterminator='\n');wr.writerow(['x','T_K','Lambda_erg_cm3_s']);wr.writerows(cool.rows)
    if not all(x['passed'] for x in checks):raise RuntimeError(json.dumps(checks,indent=2))
    with (out/'reference-comparison.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(readings[0]),lineterminator='\n');wr.writeheader();wr.writerows(readings)
    comparisons=[]
    for quantity,wave in [('T',None)]+[('amplification',w) for w in c['wavelength_pc']]:
        pp=[r for r in readings if r['quantity']==quantity and r['wavelength_pc']==wave];scored=[r['residual_dex'] for r in pp if r['residual_dex'] is not None]
        comparisons.append(dict(quantity=quantity,wavelength_pc=wave,reference_points=len(pp),scored_points=len(scored),rms_dex=float(np.sqrt(np.mean(np.array(scored)**2))),note='Descriptive sparse-point discrepancy; not a passing numerical tolerance. Points past the source-domain stop are unscored.'))
    report=dict(scope=c['scope'],checks=checks,stopping_time_Myr=float(sol.t[-1]),stop_reason='Cooling source-data boundary' if sol.status==1 else 'Configured time limit',cooling_temperature_domain_K=[cool.tmin,cool.tmax],solutions=summary,interpolation_probe=interpolation_probe,reference_comparison=comparisons,limitations=c['limitations'],environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),code_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['reproduce.py','cooling-readings.csv','reference-readings.csv']},config_sha256=hashlib.sha256((ROOT/'config.json').read_bytes()).hexdigest())
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks=len(checks),stop=report['stopping_time_Myr'],domain=report['cooling_temperature_domain_K'],solutions=summary),indent=2))


if __name__=='__main__':main()
