"""Cold-to-warm FGH/GHF prototype. Not a certified reproduction of FGH figures.

Run beside cold_branch.py, reproduce.py, config.json and config-cold.json.
See config-two-phase.json for historical inputs and declared approximations.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cold_branch import KB, recombination
from reproduce import log_coefficient, save_figure

EV=1.602176634e-12
ME=9.1093837015e-28
ROOT=Path(__file__).resolve().parent


def ramp_rate(t, threshold, slope, cutoff=np.inf):
    """Maxwellian <sigma v> for sigma=slope*(E-threshold), truncated at cutoff.

    E and kT in eV; slope in cm2/eV. Analytic integral of E*sigma*exp(-E/kT).
    """
    theta=KB*np.asarray(t)/EV
    def tail(x):
        return np.exp(-x/theta)*theta*(x*x+2*x*theta+2*theta**2-threshold*(x+theta))
    integral=tail(threshold)-(0 if np.isinf(cutoff) else tail(cutoff))
    return np.sqrt(8*EV/(np.pi*ME))*slope*np.maximum(integral,0)/theta**1.5


def rates(t, tail):
    # GHF p.178: H-ionization triangle (13.6,0), (50,1.2e-16), (300,0).
    qi=ramp_rate(t,13.6,1.2e-16/36.4,300)-ramp_rate(t,50,1.2e-16*(1/36.4+1/250),300)
    v20=np.sqrt(2*20*EV/ME)
    qip=v20*1.2e-16*(20-13.6)/36.4
    # Recover the 20-eV excitation point from rounded GHF Eqs12-13.
    qep=(3.2e-19/EV-(20-13.6)*qip)/(20-10.2)
    qeep=4.1e-17/(20*EV)
    # Only a near-threshold ramp is constrained; two *assumed* high-E tails.
    qex=ramp_rate(t,10.2,(qep/v20)/(20-10.2),20 if tail=='cut20' else np.inf)
    return qi,qex,qip,qep,qeep


def populations(n,t,c,tail='cut20'):
    n,t=np.broadcast_arrays(np.asarray(n,float),np.asarray(t,float))
    he=c['helium_nuclei_fraction']*n;nh=n-he;z=c['zeta_s_inverse']
    a,energy=recombination(t);a2,energy2=recombination(t,2)
    qi,qex,qip,qep,qeep=rates(t,tail)
    def species(ne):
        w0=a*a2*ne**2;w1=z*a2*ne;w2=z*z
        total=w0+w1+w2
        hei=he*w0/total;heii=he*w1/total;heiii=he*w2/total
        # Eliminate nHI and solve the superthermal balance as a quadratic.
        h=nh*a*ne;A=z+ne*(qi+a);D=ne*qeep
        aa=D*qip;bb=D*A+h*(qip+qep)-z*(hei+heii)*qip
        cc=z*(h+(hei+heii)*A);disc=np.sqrt(bb*bb+4*aa*cc)
        nf=2*cc/(bb+disc)
        hi=h/(A+qip*nf);hii=nh-hi
        return hi,hii,hei,heii,heiii,nf
    lo=np.zeros_like(n);hi=nh+2*he
    for _ in range(64):
        ne=(lo+hi)/2;h0,h1,he0,he1,he2,nf=species(ne)
        over=ne>h1+he1+2*he2
        hi=np.where(over,ne,hi);lo=np.where(over,lo,ne)
    ne=(lo+hi)/2;h0,h1,he0,he1,he2,nf=species(ne)
    return dict(ne=ne,HI=h0,HII=h1,HeI=he0,HeII=he1,HeIII=he2,nf=nf,NH=nh,
                alpha=a,alpha2=a2,rec_energy=energy,rec_energy2=energy2,
                qi=qi,qex=qex,qip=qip,qep=qep,qeep=qeep)


def balance(n,t,c,ab,table,tail='cut20',rec_factor=1.):
    p=populations(n,t,c,tail);ne=p['ne'];z=c['zeta_s_inverse']
    heat=p['nf']*(3.2e-19*p['HI']+4.1e-17*ne)+4.6e-10*ne*z
    metal=np.zeros_like(ne)
    for row in table['transitions']:
        element={'cii':'C','siii':'Si','feii':'Fe','oi':'O'}[row['id'].split('-')[0]]
        for collider,density in [('electron',ne),('hydrogen',p['HI'])]:
            if row[collider] is not None:
                metal+=p['NH']*ab[element]*density*np.exp(log_coefficient(t,row,collider))
    hydrogen=ne*p['HI']*EV*(13.6*p['qi']+10.2*p['qex'])
    rec=rec_factor*ne*((p['HII']+p['HeII'])*p['rec_energy']+p['HeIII']*p['rec_energy2'])
    ff=1.5e-27*np.sqrt(t)*ne*(p['HII']+p['HeII']+4*p['HeIII'])
    cool=metal+hydrogen+rec+ff
    return (cool-heat)/(cool+heat),dict(**p,heat=heat,cool=cool,metal=metal,hydrogen=hydrogen,rec=rec,ff=ff)


def curve(ts,c,ab,table,tail='cut20',rec_factor=1.,samples=None):
    # Parameterize by T to retain folds in T(n), rather than drawing across them.
    ns=np.geomspace(c['n0_min_cm-3'],c['n0_max_cm-3'],samples or c['density_scan_samples'])
    f,_=balance(ns[None,:],ts[:,None],c,ab,table,tail,rec_factor)
    ii,jj=np.where(f[:,:-1]*f[:,1:]<0)
    t=ts[ii];lo=ns[jj];hi=ns[jj+1];flo=f[ii,jj]
    for _ in range(42):
        mid=np.sqrt(lo*hi);fm,_=balance(mid,t,c,ab,table,tail,rec_factor)
        same=fm*flo>0;lo=np.where(same,mid,lo);hi=np.where(same,hi,mid);flo=np.where(same,fm,flo)
    n=np.sqrt(lo*hi);err,p=balance(n,t,c,ab,table,tail,rec_factor)
    return n,t,err,p


def main():
    c=json.loads((ROOT/'config-two-phase.json').read_text());base=json.loads((ROOT/'config-cold.json').read_text());table=json.loads((ROOT/'config.json').read_text())
    out=ROOT/'output';out.mkdir(exist_ok=True);checks=[];summaries=[];rows=[]
    def check(name,value,tol):
        checks.append(dict(name=name,error=float(value),tolerance=tol,passed=bool(value<tol)))
    def condition(name,value):checks.append(dict(name=name,passed=bool(value)))
    # Independent numerical Maxwellian quadrature of the H ionization triangle.
    e=np.linspace(13.6,300,50001);sigma=np.where(e<=50,1.2e-16*(e-13.6)/36.4,1.2e-16*(300-e)/250)
    for t in [3000.,8000.,12000.]:
        theta=KB*t/EV
        numeric=np.sqrt(8*EV/(np.pi*ME))/theta**1.5*np.trapezoid(e*sigma*np.exp(-e/theta),e)
        check('H ionization quadrature at '+str(t),abs(rates(t,'cut20')[0]/numeric-1),1e-4)
    fig,axes=plt.subplots(1,3,figsize=(15,5.6),layout='constrained')
    ns=np.geomspace(c['T_max_K'],c['T_min_K'],c['temperature_samples'])
    for (case,div),color in zip(base['depletion_divisors'].items(),['#22645c','#8b7050','#776b80']):
        ab={e:v/div[e] for e,v in base['abundances_per_H_nucleus'].items()}
        n,t,err,p=curve(ns,c,ab,table)
        condition(case+': nonempty unique root at sampled temperatures',len(t)>0 and len(np.unique(t))==len(t))
        check(case+': thermal balance',np.max(abs(err)),1e-9)
        check(case+': charge balance',np.max(abs((p['HII']+p['HeII']+2*p['HeIII'])/p['ne']-1)),1e-10)
        z=c['zeta_s_inverse'];ne=p['ne'];nf=p['nf']
        production=z*(p['HI']+p['HeI']+p['HeII'])
        removal=nf*(p['HI']*(p['qip']+p['qep'])+ne*p['qeep'])
        check(case+': superthermal balance',np.max(abs(removal/production-1)),1e-10)
        ion=p['HI']*(z+ne*p['qi']+nf*p['qip'])
        check(case+': H ionization',np.max(abs(p['alpha']*ne*p['HII']/ion-1)),1e-9)
        check(case+': He I ionization',np.max(abs(p['alpha']*ne*p['HeII']/(z*p['HeI'])-1)),1e-10)
        check(case+': He II ionization',np.max(abs(p['alpha2']*ne*p['HeIII']/(z*p['HeII'])-1)),1e-10)
        # GHF Eq21 is an independent conservation residual after eliminating Eq20.
        electron_in=nf*p['HI']*(2*p['qip']+p['qep'])+nf*ne*p['qeep']+ne*p['HI']*p['qi']
        electron_out=ne*(p['alpha']*(p['HII']+p['HeII'])+p['alpha2']*p['HeIII'])
        check(case+': thermal-electron balance (Eq21)',np.max(abs(electron_in/electron_out-1)),1e-9)
        n2,t2,_,_=curve(ns,c,ab,table,samples=2*c['density_scan_samples']-1)
        condition(case+': scan refinement root counts',np.array_equal(t,t2))
        check(case+': scan refinement densities',np.max(abs(n/n2-1)),1e-9)
        pressure=(n+ne)*t
        slope=np.diff(np.log(pressure))/np.diff(np.log(n))
        condition(case+': pressure curve has a negative-slope interval',np.any(slope<0))
        # These are turning points of this prototype, not measured FGH values.
        dp=np.diff(np.log(pressure));turns=np.flatnonzero(dp[:-1]*dp[1:]<0)
        summaries.append(dict(case=case,roots=len(n),sampled_temperatures_without_bracket=len(ns)-len(t),turning_point_grid_brackets=[dict(n0_cm3=[float(n[i]),float(n[i+2])],P_over_k=[float(pressure[i]),float(pressure[i+2])]) for i in turns]))
        for ax,y in zip(axes,[t,pressure,ne]):ax.loglog(n,y,color=color,label=case,lw=1.8)
        # Show declared model sensitivity; do not call this a statistical uncertainty band.
        for tail,recfactor,ls in [('linear',1.,':'),('cut20',2/3,'--')]:
            nv,tv,ev,pv=curve(ns,c,ab,table,tail,recfactor)
            condition(case+f': sensitivity {tail}/{recfactor} roots',np.array_equal(t,tv))
            check(case+f': sensitivity {tail}/{recfactor} thermal residual',np.max(abs(ev)),1e-9)
            summaries[-1][f'max_relative_density_change_at_fixed_T_{tail}_{recfactor}']=float(np.max(abs(nv/n-1)))
            for ax,y in zip(axes,[tv,(nv+pv['ne'])*tv,pv['ne']]):ax.loglog(nv,y,color=color,ls=ls,lw=.9,alpha=.65)
        for i in range(len(n)):
            rows.append([case,n[i],t[i],pressure[i],ne[i],nf[i],*[p[k][i] for k in ['HI','HII','HeI','HeII','HeIII','heat','cool','metal','hydrogen','rec','ff']],err[i]])
    if not all(x['passed'] for x in checks):raise RuntimeError(json.dumps(checks,indent=2))
    for ax,ylabel,title in zip(axes,['Temperature (K)',r'$P/k=(n_0+n_e)T$ (K cm$^{-3}$)',r'$n_e$ (cm$^{-3}$)'],['Cold and warm equilibrium','The two-phase pressure structure','Thermal electrons']):
        ax.set(xlabel=r'$n_0$ (H+He nuclei cm$^{-3}$)',ylabel=ylabel,title=title);ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0].legend(frameon=False)
    fig.suptitle('Toward FGH (1969): cold-to-warm equilibrium prototype')
    fig.supxlabel('Solid: declared baseline. Dotted: alternative Lyman-alpha tail. Dashed: recombination-energy factor 2/3.\nHistorical inputs with explicit approximations; not a certified reproduction of FGH Figures 1 and 2.',fontsize=9)
    matplotlib.rcParams['svg.hashsalt']='fgh-two-phase'
    save_figure(fig,out,'two-phase');plt.close(fig)
    with (out/'two-phase.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['case','n0_cm-3','T_K','P_over_k_K_cm-3','ne_cm-3','n_superthermal_cm-3','HI_cm-3','HII_cm-3','HeI_cm-3','HeII_cm-3','HeIII_cm-3','heating_erg_cm-3_s-1','cooling_erg_cm-3_s-1','metal_cooling','H_excitation_ionization_cooling','recombination_cooling','free_free_cooling','normalized_thermal_residual']);w.writerows(rows)
    report=dict(scope=c['scope'],checks=checks,solutions=summaries,limitations=c['limitations'],environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),code_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['two_phase.py','cold_branch.py','reproduce.py','config.json','config-cold.json']},config_sha256=hashlib.sha256((ROOT/'config-two-phase.json').read_bytes()).hexdigest())
    (out/'validation-two-phase.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),solutions=summaries),indent=2))


if __name__=='__main__':main()
