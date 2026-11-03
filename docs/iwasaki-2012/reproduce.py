"""Generate auditable computed prototypes toward Iwasaki & Inutsuka Figs 5/7."""
from __future__ import annotations
import argparse,csv,hashlib,json,platform,time
from collections import Counter
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from front_model import *
plt.rcParams['svg.hashsalt']='iwasaki-2012-front-prototype'


def validate(audit,config,psat):
    checks=[];cases=[]
    def check(name,error,tolerance):checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,value,**detail):checks.append(dict(name=name,passed=bool(value),**detail))
    check('Cold turning point agrees with printed 1596/1598 pressure range',min(abs(P_MIN-1596),abs(P_MIN-1598)),2.)
    check('Warm turning point near printed Pmax=5012; difference retained',abs(P_MAX-5012),6.)
    check('Saturation pressure near printed 2823 K cm^-3',abs(psat-2823),1.)
    check('Conductivity-weighted static-front integral closes',abs(saturation_integral(psat)),1e-6)
    zp={p:zp_flux(p,psat) for p in [3000,4000]}
    check('ZP flux at P=3000 agrees with rounded printed 0.0032',abs(zp[3000]/.0032-1),.02)
    check('Independent left/right phase integrations match for ZP flux',abs(zp_mismatch(3000,zp[3000])),1e-5)
    check('ZP eigenvalue stable under tighter integration tolerance',abs(zp_flux(3000,psat,rtol=7e-9)/zp[3000]-1),1e-5)
    for c in audit['classification_examples']:
        p=c['P_c_over_kB'];j=c['j_over_j_ZP']*zp[p];front=Front(p,j)
        code,detail=classify(front,zp[p],config['ode_relative_tolerance'],config['CNM_temperature_offset_fraction'],config['maximum_integration_temperature_K'])
        code2,_=classify(front,zp[p],config['ode_relative_tolerance']/4,config['CNM_temperature_offset_fraction']/2,config['maximum_integration_temperature_K'])
        condition(f"Published case P={p}, j/jZP={c['j_over_j_ZP']}: {c['expected']}",LABELS[code]==c['expected'],computed=LABELS[code],source=c['source'])
        condition(f"Refinement preserves case P={p}, ratio={c['j_over_j_ZP']}",code2==code)
        cases.append(dict(**c,minus_j_normalized=j,computed=LABELS[code],details=detail))
    front=Front(4000,zp[4000]*500);ts=np.linspace(front.tc,front.tcrit*(1-1e-5),50);n=front.density_number(ts)
    rho=MH*n;v=front.j/rho
    check('Mass flux conserved in reconstructed density/velocity',np.max(abs(rho*v/front.j-1)),1e-12)
    check('Momentum flux conserved in reconstructed states',np.max(abs((rho*v*v+front.pressure(ts))/front.momentum-1)),1e-12)
    check('Critical-point pressure is exactly half the momentum flux',abs(front.pressure(front.tcrit)/(front.momentum/2)-1),1e-7)
    h=ts*1e-5;analytic=front.cooling_derivative(ts)
    numerical=(front.cooling(ts+h)-front.cooling(ts-h))/(2*h)
    # Exclude very near the critical point where a centered step crosses it.
    check('Loss derivative agrees with independent finite differences',np.max(abs(analytic[:-2]-numerical[:-2]))/np.max(abs(analytic[:-2])),1e-6)
    condition('No outgoing saddle is assigned a physical class',classify(Front(1700,50))[0]==5)
    if not all(c['passed'] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    return checks,cases,zp


def save(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=180);p=out/(name+'.svg');fig.savefig(p,metadata={'Date':None})
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n');plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--audit',type=Path,default=Path(__file__).with_name('source-audit.json'))
    parser.add_argument('--output',type=Path,default=Path('output'))
    args=parser.parse_args();start=time.perf_counter();config=json.loads(args.config.read_text());audit=json.loads(args.audit.read_text());out=args.output;out.mkdir(parents=True,exist_ok=True)
    psat=saturation_pressure();checks,cases,zp=validate(audit,config,psat)
    print(f'Validation: {len(checks)} checks passed; Psat={psat:.6f}; jZP(3000)={zp[3000]:.9g}',flush=True)
    temperatures=np.geomspace(10,20000,1201);densities=1/cooling_ratio(temperatures);pressures=temperatures*densities
    with (out/'equilibrium.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['T_K','nH_cm-3','P_over_kB_K_cm-3']);w.writerows(zip(temperatures,densities,pressures))
    fig,ax=plt.subplots(figsize=(8.5,4.7),layout='constrained');ax.loglog(densities,pressures,color='#22645c',lw=1.8)
    for p,label in [(P_MIN,'Pmin'),(psat,'Psat'),(P_MAX,'Pmax')]:ax.axhline(p,color='#8b7050',lw=.8);ax.text(.06,p*1.035,f'{label} / kB = {p:.2f}',fontsize=8)
    ax.set(xlim=(.04,1000),ylim=(900,15000),xlabel=r'H-atom density n (cm$^{-3}$)',ylabel=r'Pressure P/kB (K cm$^{-3}$)',title='Equilibrium landmarks for the front calculation')
    ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    fig.supxlabel('Printed cooling fit, Eqs (5)–(7). Pmin/Pmax differ slightly from the rounded labels in the source.',fontsize=8)
    save(fig,out,'equilibrium')
    colors=['#c9b792','#9cb8aa','#b0a0bc','#3e7168','#293434','#e0dfd9','#ffffff']
    cmap=ListedColormap(colors);norm=BoundaryNorm(np.arange(-.5,7.5),cmap.N);grid_summary=[]
    for name in ['figure5','figure7']:
        settings=config[name]
        makep=np.linspace if name=='figure5' else np.geomspace
        ps=makep(settings['pressure_min_K_cm_inverse3'],settings['pressure_max_K_cm_inverse3'],settings['pressure_samples'])
        js=np.geomspace(settings['minus_j_min_normalized'],settings['minus_j_max_normalized'],settings['flux_samples'])
        grid=np.full((len(js),len(ps)),6,dtype=int);rows=[];boundaries=[];critical=[]
        for column,p in enumerate(ps):
            jzp=zp_flux(float(p),psat) if psat<p<P_MAX else 0
            boundaries.append(jzp if jzp>0 else np.nan)
            tc=equilibrium_temperatures(float(p))[0];nc=p/tc
            critical.append(nc*np.sqrt(R)*tc/(np.sqrt(10000)+np.sqrt(10000-tc))/1e5)
            for row,j in enumerate(js):
                front=Front(float(p),float(j));code,detail=classify(front,jzp,config['ode_relative_tolerance'],config['CNM_temperature_offset_fraction'],config['maximum_integration_temperature_K'])
                grid[row,column]=code
                rows.append([p,j,front.j,code,LABELS[code],front.tc,front.tcrit,detail.get('T_terminal_K',''),detail.get('flux_distance_exponent',''),detail['reason']])
            if column%5==0:print(f'{name}: pressure column {column+1}/{len(ps)}',flush=True)
        with (out/(name+'-grid.csv')).open('w',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(['Pc_over_kB_K_cm-3','minus_j_normalized','j_g_cm-2_s-1','class_code','class_label','Tc_K','Tcrit_K','terminal_T_K','flux_distance_exponent','decision']);w.writerows(rows)
        with (out/(name+'-boundaries.csv')).open('w',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(['Pc_over_kB_K_cm-3','minus_jZP_normalized','minus_j_at_Tcrit_10000K_normalized']);w.writerows(zip(ps,boundaries,critical))
        fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
        ax.pcolormesh(ps,js,grid,cmap=cmap,norm=norm,shading='nearest',rasterized=True)
        ax.plot(ps,boundaries,'--',color='white',lw=1.5,label='condensation ZP boundary')
        if name=='figure7':ax.plot(ps,critical,'--',color='#413d42',lw=1,label='Tcrit = 10,000 K')
        ax.axvline(psat,color='#606b6b',ls=':',lw=1)
        ax.set(yscale='log',xlabel=r'CNM pressure $P_c/k_B$ (K cm$^{-3}$)',ylabel=r'Condensation flux $-j$ / [$m_H$ cm$^{-3}$ km s$^{-1}$]',ylim=(js[0],js[-1]),xlim=(ps[0],ps[-1]),title=f'Toward Iwasaki & Inutsuka (2012), Figure {name[-1]}: computed subset')
        if name=='figure7':ax.set_xscale('log')
        handles=[Patch(facecolor=c,edgecolor='#606b6b',label=LABELS[i]) for i,c in enumerate(colors)]
        ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),frameon=False,fontsize=8)
        fig.supxlabel(f'{len(ps)} × {len(js)} computed samples. Gray: outgoing CNM branch not available. White: unresolved or beyond the T window; no interpolation of classes.\nZP boundary from matched phase integrations. The published jcrit curve is omitted pending a reference-state audit.',fontsize=8.5)
        save(fig,out,name)
        counts={LABELS[i]:int(np.sum(grid==i)) for i in range(7)};grid_summary.append(dict(figure=name,pressure_samples=len(ps),flux_samples=len(js),classification_counts=counts,decision_counts=dict(Counter(r[-1] for r in rows))))
    report=dict(scope=config['scope'],checks=checks,published_examples=cases,grid_summary=grid_summary,
        calculated_landmarks=dict(Pmin_over_kB=P_MIN,Pmax_over_kB=P_MAX,Psat_over_kB=psat,minus_jZP_at_P3000=zp[3000]),
        limitations=config['assumptions'],open_questions_for_review=audit['open_questions_for_review'],
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        elapsed_seconds=time.perf_counter()-start,code_sha256={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['front_model.py','reproduce.py','source-audit.json']},config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(checks_passed=len(checks),grids=grid_summary,elapsed_seconds=report['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
