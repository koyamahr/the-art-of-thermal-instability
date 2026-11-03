"""Compare a declared FGH prototype against 96 sparse readings of original figures.

Run beside two_phase.py and its dependencies, digitized-fgh.csv, and
config-comparison.json. No source PDF is needed: all pixel readings and axis
calibrations are included. Reading errors are envelopes, not statistical errors.
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
from two_phase import balance, curve, EV
from reproduce import log_coefficient, save_figure

ROOT=Path(__file__).resolve().parent


def calibrate(ticks):
    pixels,values=np.array(ticks,dtype=float).T
    slope,intercept=np.polyfit(pixels,values,1)
    residual=float(np.max(abs(slope*pixels+intercept-values)))
    return slope,intercept,residual


def read_points(config):
    result=[]
    for row in csv.DictReader((ROOT/'digitized-fgh.csv').open()):
        fig=config['figures'][row['figure']]
        sx,bx,ex=calibrate(fig['x_ticks']);sy,by,ey=calibrate(fig['y_ticks'])
        x=float(row['x_px']);y=float(row['y_px'])
        logn=sx*x+bx;logvalue=sy*y+by-(3 if row['quantity']=='ne' else 0)
        result.append(dict(**row,n=10**logn,value=10**logvalue,
                           logn_envelope=abs(sx)*float(row['reading_halfwidth_x_px'])+ex,
                           logvalue_envelope=abs(sy)*float(row['reading_halfwidth_y_px'])+ey))
    return result


def measured_balance(n,t,c,ab,table,options=None):
    """One-factor diagnostics are deliberately separate from the published model.

    Coefficient multipliers are sensitivity probes, not recovered rates or fits.
    In particular a thermal H excitation probe does not change the 20-eV point.
    """
    options=options or {}
    _,p=balance(n,t,c,ab,table,rec_factor=options.get('rec_factor',1.))
    indirect=p['nf']*(3.2e-19*p['HI']+4.1e-17*p['ne'])
    direct=4.6e-10*p['ne']*c['zeta_s_inverse']
    heat=indirect*options.get('indirect_scale',1.)+direct*options.get('direct_scale',1.)
    metal=np.zeros_like(p['ne']);terms={}
    for row in table['transitions']:
        el={'cii':'C','siii':'Si','feii':'Fe','oi':'O'}[row['id'].split('-')[0]]
        value=np.zeros_like(p['ne'])
        for collider,density in [('electron',p['ne']),('hydrogen',p['HI'])]:
            if row[collider] is not None:
                value+=p['NH']*ab[el]*density*np.exp(log_coefficient(t,row,collider))
        value*=options.get(row['id']+'_scale',1.)
        terms[row['id']]=value;metal+=value
    lya=p['ne']*p['HI']*10.2*EV*p['qex']*options.get('lya_scale',1.)
    hion=p['ne']*p['HI']*13.6*EV*p['qi']
    cool=metal+lya+hion+p['rec']+p['ff']
    budget=dict(indirect_heating=indirect*options.get('indirect_scale',1.),direct_heating=direct*options.get('direct_scale',1.),
                **terms,H_excitation=lya,H_ionization=hion,recombination=p['rec'],free_free=p['ff'])
    return (cool-heat)/(cool+heat),dict(**p,diagnostic_heat=heat,diagnostic_cool=cool,budget=budget)


def roots_at_temperature(ts,c,ab,table,options=None,samples=161):
    ts=np.asarray(ts,float)
    ns=np.geomspace(c['n0_min_cm-3'],c['n0_max_cm-3'],samples)
    f,_=measured_balance(ns[None,:],ts[:,None],c,ab,table,options)
    ii,jj=np.where(f[:,:-1]*f[:,1:]<0)
    t=ts[ii];lo=ns[jj];hi=ns[jj+1];flo=f[ii,jj]
    for _ in range(44):
        mid=np.sqrt(lo*hi);fm,_=measured_balance(mid,t,c,ab,table,options)
        same=fm*flo>0;lo=np.where(same,mid,lo);hi=np.where(same,hi,mid);flo=np.where(same,fm,flo)
    n=np.sqrt(lo*hi);err,p=measured_balance(n,t,c,ab,table,options)
    return ii,n,t,err,p


def roots_at_density(ns,c,ab,table,samples=241):
    ns=np.asarray(ns,float);ts=np.geomspace(c['T_min_K'],c['T_max_K'],samples)
    f,_=balance(ns[:,None],ts[None,:],c,ab,table)
    ii,jj=np.where(f[:,:-1]*f[:,1:]<0)
    n=ns[ii];lo=ts[jj];hi=ts[jj+1];flo=f[ii,jj]
    for _ in range(44):
        mid=np.sqrt(lo*hi);fm,_=balance(n,mid,c,ab,table)
        same=fm*flo>0;lo=np.where(same,mid,lo);hi=np.where(same,hi,mid);flo=np.where(same,fm,flo)
    t=np.sqrt(lo*hi);err,p=balance(n,t,c,ab,table)
    return ii,n,t,err,p


def warm_peak(c,ab,table):
    def pressure(logt):
        _,n,t,_,p=roots_at_temperature([np.exp(logt)],c,ab,table)
        if len(n)!=1:raise ValueError('Warm-peak search requires one density root')
        return float((n[0]+p['ne'][0])*t[0])
    # The source anchor is not used: bracket a model maximum from a coarse scan.
    ts=np.geomspace(5000,11000,71)
    _,n,t,_,p=roots_at_temperature(ts,c,ab,table)
    ps=(n+p['ne'])*t;i=int(np.argmax(ps))
    if i in (0,len(t)-1):raise ValueError('Warm peak is not bracketed')
    lo,hi=np.log(t[i-1]),np.log(t[i+1]);ratio=(np.sqrt(5)-1)/2
    x1=hi-ratio*(hi-lo);x2=lo+ratio*(hi-lo);f1=pressure(x1);f2=pressure(x2)
    for _ in range(32):
        if f1<f2:
            lo=x1;x1=x2;f1=f2;x2=lo+ratio*(hi-lo);f2=pressure(x2)
        else:
            hi=x2;x2=x1;f2=f1;x1=hi-ratio*(hi-lo);f1=pressure(x1)
    temp=np.exp((lo+hi)/2);_,n,t,err,p=roots_at_temperature([temp],c,ab,table)
    return dict(n0_cm3=float(n[0]),T_K=float(t[0]),P_over_k=float((n[0]+p['ne'][0])*t[0]),ne_cm3=float(p['ne'][0]),thermal_residual=float(err[0]),logT_bracket_width=float(hi-lo))


def main():
    config=json.loads((ROOT/'config-comparison.json').read_text());c=json.loads((ROOT/'config-two-phase.json').read_text())
    base=json.loads((ROOT/'config-cold.json').read_text());table=json.loads((ROOT/'config.json').read_text())
    points=read_points(config);out=ROOT/'output';out.mkdir(exist_ok=True)
    checks=[];results=[];summary=[];maxres=[]
    def check(name,value,tolerance):checks.append(dict(name=name,error=float(value),tolerance=tolerance,passed=bool(value<tolerance)))
    def condition(name,value):checks.append(dict(name=name,passed=bool(value)))
    condition('Reference set contains exactly96 original-figure readings',len(points)==96)
    for name,fig in config['figures'].items():
        for axis in ['x','y']:
            slope,offset,residual=calibrate(fig[axis+'_ticks'])
            check(name+' '+axis+': calibration residual in dex',residual,.01)
            condition(name+' '+axis+': axis orientation',slope>0 if axis=='x' else slope<0)
    fig,axes=plt.subplots(1,3,figsize=(15,5.8),layout='constrained')
    for (case,div),color in zip(base['depletion_divisors'].items(),['#22645c','#8b7050','#776b80']):
        ab={e:v/div[e] for e,v in base['abundances_per_H_nucleus'].items()}
        ref=[x for x in points if x['case']==case and x['quantity']=='T']
        condition(case+': source temperature points are monotonic in density',np.all(np.diff([x['n'] for x in ref])>0))
        ts=np.array([x['value'] for x in ref])
        ii,n,t,err,p=roots_at_temperature(ts,c,ab,table)
        condition(case+': one density root per reference temperature',np.array_equal(ii,np.arange(len(ts))))
        maxres.extend(abs(err));ii2,n2,t2,_,_=roots_at_temperature(ts,c,ab,table,samples=321)
        condition(case+': doubled density scan retains root count',np.array_equal(ii,ii2))
        check(case+': density scan refinement',np.max(abs(n/n2-1)),1e-9)
        # Verify this independent comparison wrapper leaves the original model unchanged.
        nr,tr,er,pr=curve(ts,c,ab,table)
        check(case+': baseline agrees with published prototype density roots',np.max(abs(n/nr-1)),1e-9)
        for j,x in enumerate(ref):
            delta=np.log10(n[j]/x['n'])
            results.append(dict(**x,comparison='density at reference temperature',model_value=float(n[j]),residual_dex=float(delta),candidate_values=[float(n[j])],root_count=1))
        for quantity in ['P','ne']:
            refq=[x for x in points if x['case']==case and x['quantity']==quantity]
            ns=np.array([x['n'] for x in refq]);ix,nq,tq,eq,pq=roots_at_density(ns,c,ab,table)
            ix2,nq2,tq2,_,_=roots_at_density(ns,c,ab,table,samples=481)
            condition(case+' '+quantity+': doubled temperature scan root counts',np.array_equal(ix,ix2))
            check(case+' '+quantity+': temperature scan refinement',np.max(abs(tq/tq2-1)),1e-9)
            maxres.extend(abs(eq))
            values=(nq+pq['ne'])*tq if quantity=='P' else pq['ne']
            for j,x in enumerate(refq):
                candidates=values[ix==j].tolist()
                # Do not select the closest branch to improve the error score.
                v=candidates[0] if len(candidates)==1 else None
                results.append(dict(**x,comparison='value at reference density',model_value=v,residual_dex=float(np.log10(v/x['value'])) if v is not None else None,candidate_values=candidates,root_count=len(candidates)))
        # Only draw depleted model branches over the density interval present in the source.
        nt,tt,_,pt=curve(np.geomspace(c['T_max_K'],c['T_min_K'],501),c,ab,table)
        for ax,quantity,ys in zip(axes,['T','P','ne'],[tt,(nt+pt['ne'])*tt,pt['ne']]):
            refs=[x for x in points if x['case']==case and x['quantity']==quantity]
            mask=(nt>=min(x['n'] for x in refs))&(nt<=max(x['n'] for x in refs))
            ax.loglog(np.where(mask,nt,np.nan),np.where(mask,ys,np.nan),color=color,lw=1.6,label=case)
            for x in refs:
                ax.errorbar(x['n'],x['value'],xerr=[[x['n']*(1-10**(-x['logn_envelope']))],[x['n']*(10**x['logn_envelope']-1)]],
                            yerr=[[x['value']*(1-10**(-x['logvalue_envelope']))],[x['value']*(10**x['logvalue_envelope']-1)]],fmt='o',mfc='white',mec=color,ecolor=color,ms=3.5,lw=.6,alpha=.85)
        for q in ['T','P','ne']:
            rr=[x for x in results if x['case']==case and x['quantity']==q];scored=[x['residual_dex'] for x in rr if x['residual_dex'] is not None]
            summary.append(dict(case=case,quantity=q,metric='log10(n_model(T_ref)/n_ref)' if q=='T' else 'log10(y_model(n_ref)/y_ref)',points=len(rr),scored=len(scored),unscored=len(rr)-len(scored),median_signed_dex=float(np.median(scored)),rms_dex=float(np.sqrt(np.mean(np.array(scored)**2))),max_absolute_dex=float(max(abs(np.array(scored))))))
    check('Normalized thermal residual at all scored roots',max(maxres),1e-9)
    normal=base['abundances_per_H_nucleus'];peak=warm_peak(c,normal,table)
    check('Refined warm-peak bracket in logT',peak['logT_bracket_width'],1e-7)
    anchor=config['printed_anchor'];f,p=measured_balance(anchor['n0_cm3'],anchor['T_K'],c,normal,table)
    heat=float(p['diagnostic_heat']);cool=float(p['diagnostic_cool'])
    budget={k:float(v) for k,v in p['budget'].items()}
    energy=dict(at_source_text_point=dict(n0_cm3=anchor['n0_cm3'],T_K=anchor['T_K']),ne_cm3=float(p['ne']),heating=heat,cooling=cool,cooling_over_heating=cool/heat,normalized_residual=float(f),terms=budget)
    check('Heating budget sums',abs((budget['direct_heating']+budget['indirect_heating'])/heat-1),1e-12)
    check('Cooling budget sums',abs(sum(v for k,v in budget.items() if not k.endswith('heating'))/cool-1),1e-12)
    # Check the independently printed approximate heating reduction (GHF Eq23).
    x=np.geomspace(1e-4,1,101)
    from two_phase import rates
    _,_,qip,qep,qeep=rates(7800.,'cut20')
    exact=(3.2e-19+4.1e-17*x)/(qip+qep+x*qeep)
    rounded=1.4e-11*(1+130*x)/(1+57*x)
    check('Eq13+20 heating agrees with rounded Eq23 within3percent',np.max(abs(exact/rounded-1)),.03)
    probe_specs=config['sensitivity_probes'];probes=[]
    ts=np.geomspace(5000,11000,151)
    for spec in probe_specs:
        ii,n,t,err,p=roots_at_temperature(ts,c,normal,table,spec['options'])
        condition(spec['id']+': complete warm scan',np.array_equal(ii,np.arange(len(ts))))
        check(spec['id']+': thermal residual',np.max(abs(err)),1e-9)
        pressure=(n+p['ne'])*t;j=int(np.argmax(pressure))
        condition(spec['id']+': interior warm pressure maximum',0<j<len(t)-1)
        probes.append(dict(id=spec['id'],options=spec['options'],n0_cm3=float(n[j]),T_K=float(t[j]),P_over_k=float(pressure[j]),relative_to_refined_baseline=float(pressure[j]/peak['P_over_k']-1)))
    # These are discrepancies, not tolerances that certify a historical reproduction.
    anchor_comparison=dict(model_warm_peak=peak,source_text=anchor,relative_differences={k:peak[k]/anchor[k]-1 for k in ['n0_cm3','T_K','P_over_k']})
    for ax,ylabel,title in zip(axes,['Temperature (K)',r'$P/k$ (K cm$^{-3}$)',r'$n_e$ (cm$^{-3}$)'],['FGH Figure1','FGH Figure2: pressure','FGH Figure2: electrons']):
        ax.set(xlabel=r'$n_0$ (cm$^{-3}$)',ylabel=ylabel,title=title,xlim=(.04,110));ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylim(10,15000)
    axes[1].set_ylim(50,10000)
    axes[0].legend(frameon=False)
    axes[1].plot(anchor['n0_cm3'],anchor['P_over_k'],'k*',ms=9,label='Text anchor')
    axes[1].legend(frameon=False,fontsize=8)
    fig.suptitle('FGH (1969): original-figure readings versus the declared prototype')
    fig.supxlabel('Open points:96 sparse readings with pixel-based envelopes. Lines: unchanged model; no fit to the points.\nDepleted cases are compared only where drawn in the source. Source rasters are not redistributed.',fontsize=9)
    matplotlib.rcParams['svg.hashsalt']='fgh-source-comparison';save_figure(fig,out,'comparison');plt.close(fig)
    with (out/'comparison.csv').open('w',newline='') as file:
        names=['figure','case','quantity','x_px','y_px','n','value','logn_envelope','logvalue_envelope','comparison','model_value','residual_dex','root_count','candidate_values']
        w=csv.DictWriter(file,fieldnames=names,extrasaction='ignore',lineterminator='\n');w.writeheader()
        for r in results:w.writerow({**r,'candidate_values':json.dumps(r['candidate_values'])})
    if not all(x['passed'] for x in checks):raise RuntimeError(json.dumps(checks,indent=2))
    report=dict(scope='Sparse quantitative comparison and coefficient diagnostics; historical reproduction is not certified',checks=checks,source_comparison=summary,printed_anchor_comparison=anchor_comparison,energy_budget=energy,sensitivity_probes=probes,calibration={name:{axis:dict(zip(['slope','intercept','max_residual_dex'],calibrate(spec[axis+'_ticks']))) for axis in ['x','y']} for name,spec in config['figures'].items()},method_notes=config['comparison_notes'],environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),code_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['compare.py','two_phase.py','cold_branch.py','reproduce.py','config.json','config-cold.json','config-two-phase.json','digitized-fgh.csv']},config_sha256=hashlib.sha256((ROOT/'config-comparison.json').read_bytes()).hexdigest())
    (out/'validation-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),comparisons=summary,anchor=anchor_comparison,budget=energy,probes=probes),indent=2))


if __name__=='__main__':main()
