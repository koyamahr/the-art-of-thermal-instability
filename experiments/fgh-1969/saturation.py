"""Bounded ZP saturation search on the existing FGH/GHF approximate models.

Run after two_phase.py. No extrapolation beyond its declared temperature domain.
A search without a sign bracket leaves P_sat null, not an extrapolated number.
"""
from pathlib import Path
import csv,hashlib,json
import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.optimize import brentq
from two_phase import balance,populations

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    cfg=json.loads((HERE/'config-saturation.json').read_text())
    c=json.loads((HERE/'config-two-phase.json').read_text())
    base=json.loads((HERE/'config-cold.json').read_text())
    table=json.loads((HERE/'config.json').read_text())
    out=HERE/'output';out.mkdir(exist_ok=True)
    data=list(csv.DictReader((out/'two-phase.csv').open()))
    checks=[];results=[];rows=[]
    def check(name,error,tolerance):checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def atp(t,p,ab):
        t=np.asarray(t,float);n=p/t
        for _ in range(12):
            ne=populations(n,t,c)['ne'];n=p/(t*(1+ne/n))
        _,s=balance(n,t,c,ab,table)
        return s['cool']-s['heat'],n,s
    def roots(p,ab):
        t=np.geomspace(cfg['temperature_min_K'],cfg['temperature_max_K'],257);q,*_=atp(t,p,ab)
        return [brentq(lambda z:float(atp(z,p,ab)[0]),a,b,xtol=1e-7) for a,b,x,y in zip(t[:-1],t[1:],q[:-1],q[1:]) if x*y<0]
    def integral(p,ab,count):
        rr=roots(p,ab)
        if len(rr)!=3:return None,rr,None
        t=np.linspace(rr[0],rr[-1],count);q,n,s=atp(t,p,ab);f=2500*np.sqrt(t)*q
        return float(simpson(f,x=t)/simpson(abs(f),x=t)),rr,float(np.max(abs((n+s['ne'])*t/p-1)))
    for case,div in base['depletion_divisors'].items():
        ab={e:v/div[e] for e,v in base['abundances_per_H_nucleus'].items()}
        samples=[r for r in data if r['case']==case]
        warm_boundary=max(samples,key=lambda r:float(r['T_K']))
        low=float(warm_boundary['P_over_k_K_cm-3'])*1.002
        high=max(float(r['P_over_k_K_cm-3']) for r in samples if 2000<float(r['T_K'])<10000)*.97
        scanned=[];max_pressure_error=0.;max_quad_error=0.
        if low<high:
            for p in np.linspace(low,high,cfg['pressure_samples']):
                v,rr,err=integral(p,ab,cfg['quadrature_samples'])
                if v is None:
                    rows.append([case,p,len(rr),'','','']);continue
                v2,rr2,err2=integral(p,ab,2*cfg['quadrature_samples']-1)
                max_pressure_error=max(max_pressure_error,err,err2);max_quad_error=max(max_quad_error,abs(v-v2))
                scanned.append((float(p),v2));rows.append([case,p,len(rr),v,v2,err2])
        check(case+': fixed-pressure reconstruction',max_pressure_error,1e-9)
        check(case+': quadrature refinement of normalized integral',max_quad_error,2e-4)
        brackets=[(a,b) for (a,x),(b,y) in zip(scanned[:-1],scanned[1:]) if x*y<0]
        states=[];psat=None
        if brackets:
            a,b=brackets[0];coarse=brentq(lambda p:integral(p,ab,4097)[0],a,b,xtol=1e-6)
            psat=brentq(lambda p:integral(p,ab,8193)[0],a,b,xtol=1e-6)
            check(case+': saturation pressure quadrature refinement',abs(psat/coarse-1),2e-4)
            residual,rr,err=integral(psat,ab,8193)
            check(case+': saturation integral residual',abs(residual),1e-7)
            tt=np.linspace(rr[0],rr[-1],8193);qq,_,_=atp(tt,psat,ab)
            ff=2500*np.sqrt(tt)*qq;primitive=cumulative_trapezoid(ff,x=tt,initial=0)
            check(case+': static front flux squared nonnegative',max(0.,-float(primitive.min()))/simpson(abs(ff),x=tt),1e-5)
            for label,t in zip(['cold','intermediate','warm'],rr):
                q,n,s=atp(t,psat,ab)
                check(case+': thermal closure at '+label,abs(float(q/(s['cool']+s['heat']))),1e-7)
                states.append(dict(branch=label,T_K=float(t),n0_cm3=float(n),nH_nuclei_cm3=float(.9*n),ne_cm3=float(s['ne']),P_over_kB_K_cm3=float(psat)))
        results.append(dict(case=case,P_sat_over_kB_K_cm3=psat,equilibria=states,screened_pressure_interval_K_cm3=[low,high],three_root_samples=len(scanned),normalized_integral_range=[min(x for _,x in scanned),max(x for _,x in scanned)] if scanned else None,status='conditional_root' if psat else 'no_sign_bracket_in_bounded_sample',reason='Temperature range and approximate historical inputs retained; no extrapolation or fitted rates. A positive result is model-dependent; a null result is not proof of no physical static front.'))
    if not all(c['passed'] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    with (out/'saturation-search.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['case','P_over_kB_K_cm-3','roots_in_domain','normalized_integral','refined_normalized_integral','pressure_relative_error']);w.writerows(rows)
    report=dict(scope=cfg['scope'],assumptions=cfg,checks=checks,results=results,config_sha256=sha(HERE/'config-saturation.json'),code_sha256={x:sha(HERE/x) for x in ['saturation.py','two_phase.py','cold_branch.py','reproduce.py','config-two-phase.json','config-cold.json','config.json']})
    report['code_sha256']['two-phase.csv']=sha(out/'two-phase.csv')
    (out/'validation-saturation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(results,indent=2));print(f'{len(checks)} checks passed')

if __name__=='__main__':main()
