"""A declared reduced cold-gas model built from historical FGH/GHF ingredients.

Keep reproduce.py, config.json, config-cold.json, and requirements.txt together.
Run: python cold_branch.py --output output
Numerical checks validate this implementation, not the historical full model.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reproduce import heating, log_coefficient, save_figure

KB=1.380649e-16  # erg/K


def recombination(t, z=1):
    """Excited-state number and kinetic-energy coefficients (cm3/s, erg cm3/s).

    Hjellming 1966 Eqs 52, 56-59, lambda=157890/T. Note the MINUS
    one-third exponent in gamma. For Z=2, alpha_Z=Z alpha_H(T/Z²),
    and C_Z=Z³ C_H(T/Z²). These are historical approximations.
    """
    t=np.asarray(t,dtype=float)
    if np.any(t<=0):raise ValueError('Temperature must be positive')
    lam=157890/(t/z**2);a=5.197e-14
    beta=a*np.sqrt(lam)*(.4288+.5*np.log(lam)+.469*lam**(-1/3))
    beta1=a*np.sqrt(lam)*(.7776-.4896/lam+.4796/lam**2-.057/lam**3-.405/lam**4)
    gamma=a/np.sqrt(lam)*(-.0713+.5*np.log(lam)+.640*lam**(-1/3))
    gamma1=a/np.sqrt(lam)*(.7776-.9792/lam+1.4389/lam**2-.228/lam**3-2.027/lam**4)
    return z*(beta-beta1),z**3*KB*157890*(gamma-gamma1)


def populations(n0,t,config,abundances,donors):
    """Solve charge neutrality by monotone bisection, vectorized over n0,T."""
    n0,t=np.broadcast_arrays(np.asarray(n0,dtype=float),np.asarray(t,dtype=float))
    if np.any(n0<=0):raise ValueError('Density must be positive')
    nh=(1-config['helium_nuclei_fraction'])*n0
    nhe=config['helium_nuclei_fraction']*n0
    metal=nh*sum(abundances[e] for e in ['C','Si','Fe']) if donors else np.zeros_like(nh)
    zeta=config['zeta_s_inverse'];alpha,_=recombination(t);alpha2,_=recombination(t,2)
    def species(ne):
        hii=nh*zeta/(zeta+alpha*ne)
        w0=alpha*alpha2*ne**2;w1=zeta*alpha2*ne;w2=np.full_like(ne,zeta**2)
        total=w0+w1+w2
        return hii,nhe*w0/total,nhe*w1/total,nhe*w2/total
    lo=np.zeros_like(nh);hi=nh+2*nhe+metal
    for _ in range(64):
        ne=(lo+hi)/2;hii,hei,heii,heiii=species(ne)
        positive=ne>hii+heii+2*heiii+metal
        hi=np.where(positive,ne,hi);lo=np.where(positive,lo,ne)
    ne=(lo+hi)/2;hii,hei,heii,heiii=species(ne)
    return dict(ne=ne,HI=nh-hii,HII=hii,HeI=hei,HeII=heii,HeIII=heiii,metal_e=metal,NH=nh,NHe=nhe)


def balance(n0,t,config,table,abundances,donors):
    p=populations(n0,t,config,abundances,donors)
    targets=p['HI']+p['HeI']+p['HeII']
    zeta=config['zeta_s_inverse'];x=p['ne']/p['HI']
    heat=heating(x)[2]*zeta*targets+4.6e-10*p['ne']*zeta
    metal=np.zeros_like(heat)
    for row in table['transitions']:
        element={'cii':'C','siii':'Si','feii':'Fe','oi':'O'}[row['id'].split('-')[0]]
        for channel,density in [('electron',p['ne']),('hydrogen',p['HI'])]:
            if row[channel] is not None:
                metal+=p['NH']*abundances[element]*density*np.exp(log_coefficient(t,row,channel))
    _,c1=recombination(t);_,c2=recombination(t,2)
    rec=p['ne']*((p['HII']+p['HeII'])*c1+p['HeIII']*c2)
    ff=1.5e-27*np.sqrt(t)*p['ne']*(p['HII']+p['HeII']+4*p['HeIII'])
    cool=metal+rec+ff
    return (cool-heat)/(cool+heat),dict(**p,heating=heat,cooling=cool,metal_cooling=metal,recombination_cooling=rec,free_free_cooling=ff)


def roots(densities,config,table,abundances,donors,scan=None):
    """Retain every sign-changing root found within the stated T window.

    Tangential roots may evade a sign scan. The model makes no completeness
    claim outside this window; scan refinement is tested separately.
    """
    ts=np.geomspace(config['temperature_min_K'],config['temperature_max_K'],scan or config['temperature_scan_samples'])
    residual,_=balance(densities[:,None],ts[None,:],config,table,abundances,donors)
    ii,jj=np.where((residual[:,:-1]*residual[:,1:])<0)
    ns=densities[ii];lo=ts[jj];hi=ts[jj+1];flo=residual[ii,jj]
    for _ in range(48):
        mid=np.sqrt(lo*hi);fm,_=balance(ns,mid,config,table,abundances,donors)
        same=fm*flo>0
        lo=np.where(same,mid,lo);flo=np.where(same,fm,flo);hi=np.where(same,hi,mid)
    temp=np.sqrt(lo*hi);err,p=balance(ns,temp,config,table,abundances,donors)
    return ns,temp,err,p


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config-cold.json'))
    parser.add_argument('--table',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--output',type=Path,default=Path('output'))
    args=parser.parse_args();config=json.loads(args.config.read_text());table=json.loads(args.table.read_text())
    out=args.output;out.mkdir(parents=True,exist_ok=True);checks=[]
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,result):checks.append(dict(name=name,passed=bool(result)))
    # Independent pure-H quadratic solution for number equilibrium.
    pure=dict(config,helium_nuclei_fraction=0.)
    ns=np.array([1.,100.,1000.]);ts=np.array([20.,100.,1000.]);alpha,_=recombination(ts);z=config['zeta_s_inverse']
    exact=2*z*ns/(np.sqrt(z*z+4*alpha*z*ns)+z)
    measured=populations(ns,ts,pure,{},False)['ne']
    check('Pure-H charge solution agrees with independent quadratic root',np.max(abs(measured/exact-1)),1e-12)
    a1,c1=recombination(np.geomspace(10,2000,100))
    condition('Number and energy recombination coefficients positive throughout chosen window',np.all(a1>0) and np.all(c1>0))
    check('Hydrogenic Z=2 number scaling',abs(float(recombination(400,2)[0]/recombination(100)[0])-2),1e-12)
    check('Hydrogenic Z=2 energy scaling',abs(float(recombination(400,2)[1]/recombination(100)[1])-8),1e-12)
    # Closed-form low-T leading terms of Eqs 56-59 provide an asymptotic check.
    tsmall=1e-9;lam=157890/tsmall;a,c=recombination(tsmall)
    a_asym=5.197e-14*np.sqrt(lam)*(.5*np.log(lam)-.3488)
    c_asym=5.197e-14/np.sqrt(lam)*(.5*np.log(lam)-.8489)*KB*157890
    check('Cold asymptote of number coefficient',abs(float(a/a_asym)-1),1e-6)
    check('Cold asymptote of energy coefficient',abs(float(c/c_asym)-1),1e-6)
    densities=np.geomspace(config['density_min_cm_inverse3'],config['density_max_cm_inverse3'],config['density_samples'])
    fig,axes=plt.subplots(1,2,figsize=(11,5.8),layout='constrained');rows=[];summary=[]
    for (case,divisors),color in zip(config['depletion_divisors'].items(),['#22645c','#8b7050','#776b80']):
        ab={e:v/divisors[e] for e,v in config['abundances_per_H_nucleus'].items()}
        for donors,style in [(False,'--'),(True,'-')]:
            n,t,err,p=roots(densities,config,table,ab,donors)
            label=f'{case}'+(' + metal electrons' if donors else ', H/He charge only')
            condition(label+': roots found',len(t)>0)
            check(label+': normalized thermal residual',np.max(abs(err)),1e-10)
            check(label+': H nuclei conservation',np.max(abs((p['HI']+p['HII'])/p['NH']-1)),1e-12)
            check(label+': He nuclei conservation',np.max(abs((p['HeI']+p['HeII']+p['HeIII'])/p['NHe']-1)),1e-12)
            check(label+': charge conservation',np.max(abs((p['HII']+p['HeII']+2*p['HeIII']+p['metal_e'])/p['ne']-1)),1e-12)
            ah,_=recombination(t);ahe,_=recombination(t,2)
            check(label+': H ionization balance',np.max(abs(ah*p['ne']*p['HII']/(z*p['HI'])-1)),1e-10)
            check(label+': He I ionization balance',np.max(abs(ah*p['ne']*p['HeII']/(z*p['HeI'])-1)),1e-10)
            check(label+': He II ionization balance',np.max(abs(ahe*p['ne']*p['HeIII']/(z*p['HeII'])-1)),1e-10)
            n2,t2,_,_=roots(densities,config,table,ab,donors,2*config['temperature_scan_samples']-1)
            condition(label+': refined scan finds same density/root count',np.array_equal(n,n2))
            check(label+': refined scan temperatures agree',np.max(abs(t/t2-1)),1e-10)
            pressure=(n+p['ne'])*t
            axes[0].loglog(n,t,color=color,ls=style,lw=1.7,label=label)
            axes[1].loglog(n,pressure,color=color,ls=style,lw=1.7)
            for i in range(len(n)):
                rows.append([case,int(donors),n[i],t[i],p['ne'][i],pressure[i],p['HI'][i],p['HII'][i],p['HeI'][i],p['HeII'][i],p['HeIII'][i],p['metal_e'][i],p['heating'][i],p['cooling'][i],p['metal_cooling'][i],p['recombination_cooling'][i],p['free_free_cooling'][i],err[i]])
            at100=np.isclose(n,100)
            summary.append(dict(case=case,metal_donor_electrons=donors,roots=len(t),densities_without_bracket=int(len(densities)-len(np.unique(n))),temperature_at_n0_100_K=t[at100].tolist()))
    if not all(c['passed'] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    axes[0].set(ylabel='Equilibrium temperature T (K)',ylim=(10,2000),title='Cold-gas solutions within a restricted window')
    axes[1].set(ylabel=r'Pressure proxy $P/k=(n_0+n_e)T$ (K cm$^{-3}$)',title='Pressure implied by the reduced model')
    for ax in axes:
        ax.set(xlabel=r'H + He nuclei density $n_0$ (cm$^{-3}$)',xlim=(1,1000));ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0].legend(frameon=False,fontsize=7.5,loc='upper right')
    fig.suptitle('Toward FGH (1969): a reduced cold-branch experiment',fontsize=14)
    fig.supxlabel('Primary ionization only; excited-state H/He recombination; low-density metal cooling. Solid: C+Si+Fe donor electrons.\nDashed: H/He charge only. Missing curve segments have no sign-changing root in 10–2,000 K. Not an original-figure reproduction.',fontsize=8.5)
    save_figure(fig,out,'cold-branch');plt.close(fig)
    with (out/'cold-branch.csv').open('w',newline='') as file:
        w=csv.writer(file,lineterminator='\n')
        w.writerow(['case','include_C_Si_Fe_electrons','n0_H_He_nuclei_cm-3','T_K','ne_cm-3','P_over_k_K_cm-3','nHI_cm-3','nHII_cm-3','nHeI_cm-3','nHeII_cm-3','nHeIII_cm-3','n_metal_e_cm-3','heating_erg_cm-3_s-1','cooling_erg_cm-3_s-1','metal_cooling_erg_cm-3_s-1','recombination_cooling_erg_cm-3_s-1','free_free_cooling_erg_cm-3_s-1','normalized_thermal_residual'])
        w.writerows(rows)
    hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['cold_branch.py','reproduce.py','config.json']}
    report=dict(scope=config['scope'],checks=checks,solutions=summary,limitations=config['assumptions'],
                root_search='Logarithmic sign-change scan, followed by bisection. Tangential or out-of-window roots are not certified.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),code_sha256=hashes,config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/'validation-cold.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),solutions=summary),indent=2))


if __name__=='__main__':main()
