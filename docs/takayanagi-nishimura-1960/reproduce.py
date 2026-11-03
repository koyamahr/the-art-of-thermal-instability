"""Takayanagi & Nishimura (1960): printed rates and a truncated para-H2 ladder.

This calculation closes the ladder at the highest connected level listed at
that temperature in Table 1. Missing rates are unavailable, not measured zeros.
Table 2 comparisons are diagnostics, separate from numerical validation.
Run: python reproduce.py --output output
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
plt.rcParams['svg.hashsalt']='takayanagi-nishimura-1960'


def transition_data(row,config):
    """Return upper J, upward/downward q, spontaneous A, and photon energy."""
    available=[]
    for q in row['upward_cm3_s']:
        if q is None:break
        if q<=0:raise ValueError('A tabulated rate must be positive')
        available.append(q)
    if not available:raise ValueError('No connected ladder available')
    t=row['T_K']
    if t<=0:raise ValueError('Temperature must be positive')
    j=np.arange(2,2*len(available)+1,2);up=np.array(available)
    down=up*(2*j-3)/(2*j+1)*np.exp(config['deexcitation_temperature_prefactor_K']*(2*j-1)/t)
    A=config['radiative_prefactor_s_inverse']*j*(j-1)*(2*j-1)**4/(2*j+1)
    energy=config['energy_prefactor_erg']*(2*j-1)
    return j,up,down,A,energy


def populations(nH,row,config):
    """Zero-flux solution of a finite birth/death ladder, Eqs (7)-(9)."""
    if nH<0:raise ValueError('H-atom density cannot be negative')
    j,up,down,A,energy=transition_data(row,config)
    relative=np.r_[1.,np.cumprod(nH*up/(nH*down+A))]
    p=relative/relative.sum()
    return p,float(np.dot(p[1:],A*energy))


def matrix_solution(nH,row,config):
    """Independent stationary rate-matrix solve with a normalization row."""
    j,up,down,A,energy=transition_data(row,config);m=np.zeros((len(j)+1,len(j)+1))
    for i in range(len(j)):
        upward=nH*up[i];downward=nH*down[i]+A[i]
        m[i+1,i]+=upward;m[i,i]-=upward
        m[i,i+1]+=downward;m[i+1,i+1]-=downward
    rhs=np.zeros(len(j)+1);m[-1]=1;rhs[-1]=1
    return np.linalg.solve(m,rhs)


def validate(config):
    checks=[]
    def check(name,error,tolerance):checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,result):checks.append(dict(name=name,passed=bool(result)))
    rows=config['table1']['rows'];condition('34 positive tabulated upward coefficients, blanks retained as null',sum(q is not None for r in rows for q in r['upward_cm3_s'])==34 and all(q>0 for r in rows for q in r['upward_cm3_s'] if q is not None))
    _,_,_,A,energy=transition_data(rows[0],config)
    check('Equation (7), J=2 spontaneous decay rate in s^-1',abs(A[0]/2.43648e-11-1),1e-12)
    check('J=2 to 0 photon energy in erg',abs(energy[0]/7.0683e-14-1),1e-12)
    norm=[];mins=[];matrix=[];stationary=[];energy_balance=[];dilute=[];lte=[]
    for row in rows:
        j,up,down,A,e=transition_data(row,config)
        for n in config['density_grid_cm_inverse3']:
            p,cool=populations(n,row,config)
            norm.append(abs(p.sum()-1));mins.append(np.min(p))
            matrix.append(np.max(abs(p-matrix_solution(n,row,config))))
            lhs=n*up*p[:-1];rhs=(n*down+A)*p[1:]
            stationary.append(np.max(abs(lhs-rhs)/np.maximum(lhs+rhs,1e-300)))
            collisional=np.dot(n*(p[:-1]*up-p[1:]*down),e)
            energy_balance.append(abs(collisional/cool-1))
        _,cool=populations(1e-12,row,config)
        dilute.append(abs(cool/(1e-12*up[0]*e[0])-1))
        levels=np.r_[0,j];weights=(2*levels+1)*np.exp(-.5*config['deexcitation_temperature_prefactor_K']*levels*(levels+1)/row['T_K'])
        thermal=weights/weights.sum();p,_=populations(1e14,row,config)
        lte.append(np.max(abs(p-thermal)))
    check('Population normalization',max(norm),1e-12)
    condition('All level populations nonnegative',min(mins)>=0)
    check('Recursion agrees with independent constrained rate-matrix solution',max(matrix),1e-12)
    check('Stationary transition fluxes balance',max(stationary),1e-12)
    check('Net collisional energy loss equals emitted radiative energy',max(energy_balance),1e-10)
    check('Dilute limit equals nH*q_02*DeltaE_20 per molecule',max(dilute),1e-8)
    check('Dense limit approaches Boltzmann populations of the same finite ladder',max(lte),1e-6)
    p,cool=populations(0,rows[0],config)
    condition('Zero collider density leaves all molecules in J=0 with no cooling',p[0]==1 and cool==0 and np.all(p[1:]==0))
    try:populations(-1,rows[0],config);rejected=False
    except ValueError:rejected=True
    condition('Negative collider density rejected',rejected)
    if not all(c['passed'] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    return checks


def save_figure(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=180)
    path=out/(name+'.svg');fig.savefig(path,metadata={'Date':None})
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--audit',type=Path,default=Path(__file__).with_name('source-audit.json'))
    parser.add_argument('--output',type=Path,default=Path('output'))
    args=parser.parse_args();config=json.loads(args.config.read_text());audit=json.loads(args.audit.read_text())
    checks=validate(config);out=args.output;out.mkdir(parents=True,exist_ok=True)
    rows=sorted(config['table1']['rows'],key=lambda r:r['T_K']);rates=[];results=[];comparisons=[]
    for row in rows:
        j,up,down,A,e=transition_data(row,config)
        rates.extend([row['T_K'],int(J-2),int(J),u,d,a,energy] for J,u,d,a,energy in zip(j,up,down,A,e))
        for n in config['density_grid_cm_inverse3']:
            p,cooling=populations(n,row,config)
            results.append([row['T_K'],n,int(j[-1]),cooling,*p,*([None]*(5-len(p)))])
    byT={r['T_K']:r for r in rows}
    for reference in audit['table2']['rows']:
        for n,printed in zip(audit['table2']['H_atom_densities_cm_inverse3'],reference['cooling_erg_s']):
            t=reference['T_K'];row=byT[t];_,value=populations(n,row,config)
            comparisons.append([t,n,2*sum(q is not None for q in row['upward_cm3_s']),printed,value,value/printed-1])
    for name,header,data in [
        ('rotational-rates.csv',['T_K','lower_J','upper_J','q_up_cm3_s','q_down_cm3_s','A_down_s-1','DeltaE_erg'],rates),
        ('cooling.csv',['T_K','nH_atoms_cm-3','maximum_J','cooling_per_paraH2_erg_s','p_J0','p_J2','p_J4','p_J6','p_J8'],results),
        ('table2-comparison.csv',['T_K','nH_atoms_cm-3','maximum_J','Table2_printed_erg_s','truncated_model_erg_s','fractional_difference'],comparisons)]:
        with (out/name).open('w',newline='') as file:
            w=csv.writer(file,lineterminator='\n');w.writerow(header);w.writerows(data)
    colors=['#22645c','#8b7050','#776b80','#4f7187']
    fig,axes=plt.subplots(1,2,figsize=(11,5.7),layout='constrained')
    for i,color in enumerate(colors):
        selected=[r for r in rates if r[1]==2*i]
        axes[0].loglog([r[0] for r in selected],[r[3] for r in selected],'-o',ms=3,lw=1.5,color=color,label=f'J = {2*i} → {2*i+2}')
    for n,color in zip(config['density_grid_cm_inverse3'],colors):
        selected=[r for r in results if r[1]==n];ref=sorted([r for r in comparisons if r[1]==n])
        axes[1].loglog([r[0] for r in selected],[r[3] for r in selected],'-',color=color,lw=1.6,label=f'n(H) = {n:g} cm⁻³')
        axes[1].loglog([r[0] for r in ref],[r[3] for r in ref],'x',color=color,ms=4)
    axes[0].set(title='Table 1: upward H + para-H₂ excitation',ylabel=r'Rate coefficient $q_{J\to J+2}$ (cm$^3$ s$^{-1}$)',ylim=(1e-23,3e-9))
    axes[1].set(title='Cooling from the available finite level ladder',ylabel=r'Cooling per para-H$_2$ molecule (erg s$^{-1}$)',ylim=(1e-38,1e-20))
    for ax in axes:
        ax.set(xlabel='Gas temperature T (K)',xlim=(20,5000));ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('Takayanagi & Nishimura (1960): from tabulated rates to cooling',fontsize=14)
    fig.supxlabel('Rates are plotted only where Table 1 supplies values. Lines connect samples; no interpolated coefficients are used.\nCooling: lines = truncated model, crosses = printed Table 2. The maximum included J changes with available entries.',fontsize=8.5)
    save_figure(fig,out,'rates-and-cooling')
    fig,ax=plt.subplots(figsize=(9,4.8),layout='constrained')
    for n,color in zip(config['density_grid_cm_inverse3'],colors):
        selected=sorted([r for r in comparisons if r[1]==n]);ax.semilogx([r[0] for r in selected],[100*r[5] for r in selected],'-o',color=color,ms=3,lw=1.4,label=f'n(H) = {n:g} cm⁻³')
    ax.axhline(0,color='#606b6b',lw=.8);ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    ax.set(xlabel='Gas temperature T (K)',ylabel='Difference from printed Table 2 (%)',title='External comparison: where does the truncated model differ?',xlim=(20,5000),ylim=(-25,3))
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    fig.supxlabel('100 × (truncated calculation / printed cooling − 1). These differences are measurements, not passed validation tests.\nThe ambiguous 150/155 K source row is excluded. No missing rates or corrections are invented to force agreement.',fontsize=8.5)
    save_figure(fig,out,'table2-differences')
    worst=max(comparisons,key=lambda r:abs(r[5]));comparison=dict(points=len(comparisons),maximum_absolute_fractional_difference=abs(worst[5]),worst_case=dict(T_K=worst[0],nH_atoms_cm_inverse3=worst[1],maximum_J=worst[2],printed_cooling_erg_s=worst[3],calculated_cooling_erg_s=worst[4]),
        interpretation='A finite-ladder comparison, not a validated reproduction of every Table 2 entry. The largest discrepancy is retained openly, without fitting or replacing missing coefficients.')
    report=dict(scope=config['scope'],checks=checks,external_reference_comparison=comparison,limitations=config['assumptions']+audit['caveats'],
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        code_sha256={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['reproduce.py','source-audit.json']},config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),external_reference_comparison=comparison),indent=2))


if __name__=='__main__':main()
