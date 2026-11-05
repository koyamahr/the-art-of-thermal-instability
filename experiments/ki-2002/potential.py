"""Pressure-dependent static-front potentials, with independent quadrature checks."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np
from scipy.integrate import quad,simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from equilibrium import HEATING,ratio,slope,loss,roots,static_integral,saturation_pressure
HERE=Path(__file__).resolve().parent
plt.rcParams['svg.hashsalt']='ki-static-potential'

def kappa(t):return 2500*np.sqrt(t)
def theta(t):return (5000/3)*np.asarray(t)**1.5

def potential(t,p,cold):
    """V(theta(T),P), with V(theta(Tcold),P)=0, in erg^2 cm^-4 s^-2."""
    return quad(lambda z:kappa(z)*loss(z,p/z)/1e-18,cold,t,epsabs=1e-10,epsrel=1e-10,limit=150)[0]*1e-18

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=HERE/'output');a=ap.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((HERE/'config-potential.json').read_text());ps=saturation_pressure();rs=roots(ps);scale=potential(rs[1],ps,rs[0]);coordinate_scale=theta(rs[2]);checks=[];records=[];rows=[]
    def check(name,error,tolerance):checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,ok):checks.append(dict(name=name,passed=bool(ok)))
    fig,axs=plt.subplots(2,3,figsize=(13.5,8),layout='constrained',sharey=True)
    mobile,maxs=plt.subplots(3,1,figsize=(5.5,11),layout='constrained',sharey=True)
    colors=['#356b91','#b36328','#217667']
    for col,factor in enumerate(cfg['pressure_factors']):
        p=ps*factor;r=roots(p);cold,unstable,warm=r
        temperatures=np.unique(np.r_[np.geomspace(cfg['cold_temperature_margin']*cold,cfg['warm_temperature_margin']*warm,cfg['temperature_samples']),r])
        values=np.array([potential(t,p,cold) for t in temperatures]);vr=np.array([potential(t,p,cold) for t in r]);tilt=vr[-1]-vr[0]
        # Independent composite Simpson integration, distinct from adaptive quadrature.
        fine=np.linspace(cold,warm,32001);integrand=kappa(fine)*loss(fine,p/fine)
        check(f'{factor:g} Psat: independent signed-integral cancellation',abs(simpson(integrand,x=fine)-tilt)/simpson(abs(integrand),x=fine),1e-8)
        check(f'{factor:g} Psat: Simpson mesh refinement 16001 to 32001',abs(simpson(integrand[::2],x=fine[::2])-simpson(integrand,x=fine))/simpson(abs(integrand),x=fine),1e-8)
        check(f'{factor:g} Psat: potential tilt equals existing ZP integral',abs(tilt-static_integral(p)*1e-18)/scale,1e-9)
        for branch,t,v in zip(['cold','intermediate','warm'],r,vr):
            h=t*1e-4
            numeric=(potential(t+h,p,cold)-2*v+potential(t-h,p,cold))/(h*h*kappa(t)**2)
            analytic=HEATING*(p/t)/t*(slope(t)-1)/kappa(t)
            check(f'{factor:g} Psat: {branch} curvature in theta coordinate',abs(numeric/analytic-1),2e-5)
            condition(f'{factor:g} Psat: {branch} is a '+('maximum' if branch=='intermediate' else 'minimum'),(potential(t+h,p,cold)-v)*(1 if branch!='intermediate' else -1)>0 and (potential(t-h,p,cold)-v)*(1 if branch!='intermediate' else -1)>0)
        for t in [150.,800.,3000.,6000.]:
            h=t*1e-5;numeric=(potential(t+h,p,cold)-potential(t-h,p,cold))/(theta(t+h)-theta(t-h))
            check(f'{factor:g} Psat: dV/dtheta=Q at {t:g} K',abs(numeric/loss(t,p/t)-1),2e-7)
        condition(f'{factor:g} Psat: stable minimum selection',tilt<0 if factor<1 else (tilt>0 if factor>1 else abs(tilt)/scale<1e-9))
        record=dict(pressure_factor=factor,P_over_kB_K_cm3=p,equilibria=[dict(branch=b,T_K=t,n_H_cm3=p/t,theta_erg_cm1_s1=float(theta(t)),V_erg2_cm4_s2=float(v),V_normalized=float(v/scale)) for b,t,v in zip(['cold','intermediate','warm'],r,vr)],warm_minus_cold_V_erg2_cm4_s2=float(tilt),warm_minus_cold_V_normalized=float(tilt/scale),deeper_minimum='equal within numerical tolerance' if factor==1 else ('warm' if tilt<0 else 'cold'))
        records.append(record)
        for t,v in zip(temperatures,values):rows.append([factor,p,t,p/t,theta(t),v,v/scale])
        for row,x in enumerate([temperatures,theta(temperatures)/coordinate_scale]):
            ax=axs[row,col];ax.plot(x,values/scale,color='#22645c',lw=2)
            xr=r if row==0 else theta(np.array(r))/coordinate_scale
            for xx,v,color,branch in zip(xr,vr,colors,['Cold minimum','Intermediate maximum','Warm minimum']):ax.scatter(xx,v/scale,color=color,s=34,zorder=3,label=branch)
            ax.axhline(0,color='#777',ls=':',lw=.7);ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
            ax.set_ylim(-.75,2.1)
            if row==0:ax.set_xscale('log');ax.set_xlabel('Temperature T (K)');ax.set_title(f'P = {factor:g} Psat\nP/kB = {p:.3f} K cm$^{{-3}}$')
            else:ax.set_xlabel(r'$\theta(T)\,/\,\theta(T_{warm,sat})$')
        axs[0,col].text(.04,.94,f'Warm - cold = {tilt/scale:+.6f}',transform=axs[0,col].transAxes,va='top',fontsize=9)
        ax=maxs[col];ax.semilogx(temperatures,values/scale,color='#22645c',lw=2)
        ax.scatter(r,vr/scale,c=colors,s=35,zorder=3);ax.axhline(0,color='#777',ls=':',lw=.7);ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
        ax.set(xlabel='Temperature T (K)',ylabel='V / saturation barrier',ylim=(-.75,2.1),title=f'P = {factor:g} Psat; P/kB = {p:.3f} K cm$^{{-3}}$')
        ax.text(.04,.94,f'Warm - cold = {tilt/scale:+.6f}',transform=ax.transAxes,va='top',fontsize=9)
    axs[0,0].set_ylabel(r'$V / V_{barrier,sat}$');axs[1,0].set_ylabel(r'$V / V_{barrier,sat}$');axs[0,1].legend(fontsize=8,loc='lower right')
    fig.suptitle('Static-front potential: pressure selects the deeper thermal minimum',fontsize=15)
    fig.supxlabel('Top: logarithmic temperature reveals both wells. Bottom: the natural conduction coordinate reveals their strong asymmetry.\nEach cold minimum defines V = 0; all panels use one fixed saturation-barrier scale.',fontsize=10)
    for suffix in ['png','svg']:fig.savefig(out/f'potential.{suffix}',dpi=180,metadata={'Date':None} if suffix=='svg' else None)
    mobile.suptitle('Pressure-dependent thermal potential',fontsize=13)
    mobile.supxlabel('Blue: cold minimum; orange: intermediate maximum; green: warm minimum.\nCold V = 0 separately at each pressure; one shared vertical scale.',fontsize=8)
    for suffix in ['png','svg']:mobile.savefig(out/f'potential-mobile.{suffix}',dpi=180,metadata={'Date':None} if suffix=='svg' else None)
    plt.close(mobile)
    plt.close(fig)
    with (out/'potential.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['pressure_factor','P_over_kB_K_cm3','T_K','n_H_cm3','theta_erg_cm1_s1','V_erg2_cm4_s2','V_normalized']);w.writerows(rows)
    report=dict(scope=cfg['scope'],definition='theta=(5000/3) T^(3/2); V=integral_Tcold^T kappa Q dT; dV/dtheta=Q; Vcold=0 separately at each pressure',normalization=dict(V_barrier_at_saturation_erg2_cm4_s2=scale,theta_warm_at_saturation_erg_cm1_s1=float(coordinate_scale)),saturation_pressure_over_kB_K_cm3=ps,pressures=records,checks=checks,config_sha256=hashlib.sha256((HERE/'config-potential.json').read_bytes()).hexdigest(),code_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['equilibrium.py','potential.py']})
    (out/'validation-potential.json').write_text(json.dumps(report,indent=2)+'\n')
    assert all(x['passed'] for x in checks),json.dumps([x for x in checks if not x['passed']],indent=2)
    print(f'{len(checks)} potential checks passed; normalized well tilts: '+str([r['warm_minus_cold_V_normalized'] for r in records]))
if __name__=='__main__':main()
