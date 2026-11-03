"""Generate corrected KI equilibrium curves, stability labels and ZP landmarks."""
from pathlib import Path
import argparse, csv, hashlib, importlib.util, json, platform, sys
import numpy as np
import scipy
from scipy.integrate import simpson, solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from equilibrium import HEATING, KB, ratio, slope, loss, pressure, turns, roots, static_integral, saturation_pressure

HERE = Path(__file__).resolve().parent
plt.rcParams['svg.hashsalt'] = 'ki-2002-equilibrium'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write_csv(path, header, rows):
    with path.open('w', newline='') as f:
        w = csv.writer(f, lineterminator='\n'); w.writerow(header); w.writerows(rows)


def save(fig, out, name):
    fig.savefig(out/(name+'.png'), dpi=180)
    p = out/(name+'.svg'); fig.savefig(p, metadata={'Date':None})
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('output'))
    parser.add_argument('--front-model', type=Path, help='Optional live cross-check against the Iwasaki exhibit')
    args = parser.parse_args(); out = args.output; out.mkdir(parents=True, exist_ok=True)
    c = json.loads((HERE/'config.json').read_text())
    audit = json.loads((HERE/'source-audit.json').read_text())
    t = np.geomspace(c['temperature_min_K'], c['temperature_max_K'], c['temperature_samples'])
    tc, tw = turns(); t = np.unique(np.r_[t, tc, tw])
    n = 1/ratio(t); p = t*n; beta = slope(t)
    psat = saturation_pressure(); three = roots(psat)
    checks = []
    def check(name, error, tolerance):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance, passed=bool(error<tolerance)))
    def condition(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
    check('Equilibrium energy closure over temperature grid', np.max(abs(n*ratio(t)-1)), 1e-14)
    condition('Cooling is strictly increasing in displayed interval', np.all(np.diff(ratio(t))>0))
    h = 1e-5
    numeric = (np.log(ratio(t*np.exp(h)))-np.log(ratio(t*np.exp(-h))))/(2*h)
    check('Analytic cooling logarithmic slope vs centered differences', np.max(abs(numeric-beta)), 1e-7)
    check('Both turning points satisfy beta=1', max(abs(slope(tc)-1), abs(slope(tw)-1)), 1e-12)
    condition('Cold and warm branches stable, intermediate branch unstable', slope(50)>1 and slope(1000)<1 and slope(8000)>1)
    # Check dQ/dT at fixed pressure directly, not by the slope formula.
    for ts in [50., 1000., 8000.]:
        pp = pressure(ts); dt = ts*1e-5
        direct = (loss(ts+dt, pp/(ts+dt))-loss(ts-dt, pp/(ts-dt)))/(2*dt)
        expected = HEATING*(pp/ts)/ts*(slope(ts)-1)
        check(f'Fixed-pressure loss derivative at T={ts:g}', abs(direct/expected-1), 1e-7)
    # Independent implicit inversion of the equilibrium curve at fixed densities.
    for ns in [.1, 1., 10., 100.]:
        ts = brentq(lambda x: ns*HEATING*ratio(x)-HEATING, 10, 10000, xtol=1e-10)
        check(f'Independent equilibrium root at n={ns:g}', abs(ns*ratio(ts)-1), 1e-10)
    check('Static-front integral at saturation pressure', abs(static_integral(psat)), 1e-7)
    # Independent composite quadrature of a signed integral with nonzero lobes.
    grid = np.linspace(three[0], three[-1], 16001)
    integrand = 2500*np.sqrt(grid)*loss(grid, psat/grid)/1e-18
    check('Independent Simpson cancellation of conductive integral', abs(simpson(integrand,x=grid))/simpson(abs(integrand),x=grid), 1e-7)
    condition('Saturation pressure lies inside two-phase range', pressure(tc)<psat<pressure(tw))
    condition('Static integral brackets saturation with opposite signs', static_integral(psat*.99)*static_integral(psat*1.01)<0)
    for i, ts in enumerate(three):
        check(f'Saturation equilibrium root {i+1} closes thermal balance', abs(psat/ts*ratio(ts)-1), 1e-10)
    # Fixed-pressure parcel dT/dt = -Q/(n cp), cp=(5/2)k_B per atom.
    # Compare finite time responses for cold/unstable/warm equilibria.
    for label, ts in zip(['cold','unstable','warm'],three):
        n0 = psat/ts; timescale = 2.5*KB*ts/(HEATING*abs(slope(ts)-1))
        sol = solve_ivp(lambda time,y: [-loss(y[0],psat/y[0])/(2.5*KB*(psat/y[0]))],
                        (0,.1*timescale), [ts*(1+1e-5)], rtol=1e-11, atol=1e-10)
        amplification = abs((sol.y[0,-1]-ts)/(ts*1e-5))
        check(f'Isobaric parcel response on {label} branch', abs(amplification-np.exp(.1 if label=='unstable' else -.1)), 5e-5)
    refs = audit['iwasaki']['reference_values']; rt = np.array([r['T_K'] for r in refs])
    check('Frozen Iwasaki implementation cooling-law agreement', np.max(abs(ratio(rt)/[r['Lambda_over_Gamma_cm3'] for r in refs]-1)), 1e-14)
    live = None
    if args.front_model:
        spec = importlib.util.spec_from_file_location('iwasaki_front_reference', args.front_model)
        m = importlib.util.module_from_spec(spec); sys.modules[spec.name] = m; spec.loader.exec_module(m)
        live = sha(args.front_model)
        condition('Iwasaki source matches audited version', live==audit['iwasaki']['reference_code_sha256'])
        check('Live Iwasaki cooling law on full grid', np.max(abs(ratio(t)/m.cooling_ratio(t)-1)), 1e-14)
        check('Live Iwasaki loss normalization on nonequilibrium states', np.max(abs(loss(t,2*n)/m.loss(t,2*n)-1)), 1e-14)
        check('Independent Iwasaki saturation pressure agreement', abs(psat/m.saturation_pressure()-1), 1e-10)
    # Source pressure labels are rounded/internally slightly different. Record residuals,
    # rather than using those approximate labels as coefficient fitting targets.
    landmarks = [dict(label='Pmin',T_K=float(tc),n_H_cm3=float(1/ratio(tc)),P_over_kB_K_cm3=float(pressure(tc))),
                 dict(label='Pmax',T_K=float(tw),n_H_cm3=float(1/ratio(tw)),P_over_kB_K_cm3=float(pressure(tw)))]
    saturation = [dict(branch=label,T_K=float(ts),n_H_cm3=float(psat/ts),P_over_kB_K_cm3=float(psat),isobaric_stable=bool(slope(ts)>1)) for label,ts in zip(['cold','unstable','warm'],three)]
    comparison = dict(Pmin_minus_text=float(pressure(tc)-1596),Pmin_minus_figure1=float(pressure(tc)-1598),Pmax_minus_printed=float(pressure(tw)-5012),Psat_minus_printed=float(psat-2823))
    if not all(x['passed'] for x in checks):
        raise RuntimeError(json.dumps(checks,indent=2))
    labels = np.where(abs(beta-1)<1e-12,'marginal',np.where(beta>1,'stable','unstable'))
    write_csv(out/'equilibrium.csv',['T_K','n_H_cm-3','P_over_kB_K_cm-3','Lambda_erg_cm3_s-1','beta','isobaric_branch'],zip(t,n,p,HEATING*ratio(t),beta,labels))
    write_csv(out/'saturation.csv',['branch','T_K','n_H_cm-3','P_over_kB_K_cm-3','isobaric_stable'],[list(x.values()) for x in saturation])
    write_csv(out/'coefficient-comparison.csv',['T_K','corrected_Lambda_over_Gamma_cm3','arxiv_v1_printed_Lambda_over_Gamma_cm3'],zip(t,ratio(t),ratio(t,True)))
    fig, axs = plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    curves = [(t<=tc,'#356b91','Cold / stable','-'),((t>=tc)&(t<=tw),'#b36328','Isobaric unstable','--'),(t>=tw,'#217667','Warm / stable','-')]
    for mask,col,label,ls in curves:
        for ax,x,y in zip(axs,[n,n,t],[t,p,p]):ax.loglog(x[mask],y[mask],color=col,ls=ls,lw=2,label=label)
    for ax,x,y in zip(axs,[[x['n_H_cm3'] for x in saturation]]*2+[[x['T_K'] for x in saturation]],[[x['T_K'] for x in saturation]]+[[psat]*3]*2):
        ax.scatter(x,y,s=35,c=['#356b91','#b36328','#217667'],edgecolors='white',zorder=4)
    for ax in axs[1:]:ax.axhline(psat,color='#777777',lw=.8,ls=':',label='ZP saturation')
    for ax in axs:ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    axs[0].set(xlabel=r'n = rho / mH (cm$^{-3}$)',ylabel='Temperature (K)',title='Temperature vs density')
    axs[1].set(xlabel=r'n = rho / mH (cm$^{-3}$)',ylabel=r'P / kB (K cm$^{-3}$)',title='Pressure vs density',ylim=(30,3e5))
    axs[2].set(xlabel='Temperature (K)',ylabel=r'P / kB (K cm$^{-3}$)',title='Pressure vs temperature',ylim=(30,3e5))
    axs[0].legend(fontsize=8,loc='lower left')
    fig.suptitle('Koyama & Inutsuka (2002): corrected analytic fit',fontsize=14)
    fig.supxlabel('Points: cold, unstable and warm equilibria at the ZP saturation pressure. Branch labels omit finite-wavelength conduction.',fontsize=9)
    save(fig,out,'equilibrium')
    fig,axs = plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    axs[0].loglog(t,ratio(t),label='Corrected / Iwasaki 2012',color='#217667')
    axs[0].loglog(t,ratio(t,True),'--',label='Uncorrected arXiv v1',color='#b36328')
    axs[0].set(xlabel='Temperature (K)',ylabel=r'Lambda / Gamma (cm$^3$)',title='Two coefficient corrections');axs[0].legend(fontsize=8)
    axs[1].semilogx(t,beta,color='#217667');axs[1].axhline(1,color='#777',ls=':')
    axs[1].axvspan(tc,tw,color='#b36328',alpha=.12)
    axs[1].set(xlabel='Temperature (K)',ylabel=r'beta = d ln Lambda / d ln T',title='Isobaric instability: beta < 1')
    for ax in axs:ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    save(fig,out,'coefficients')
    report = dict(scope=c['scope'],config_sha256=sha(HERE/'config.json'),code_sha256={name:sha(HERE/name) for name in ['equilibrium.py','reproduce.py','source-audit.json']},
                  versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),checks=checks,
                  turning_points=landmarks,saturation_equilibria=saturation,saturation_conductivity='kappa=2500 sqrt(T) erg cm^-1 K^-1 s^-1; planar static zero-flux end states',
                  source_pressure_residuals_K_cm3=comparison,live_front_model_sha256=live)
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'{len(checks)} checks passed; P_sat/kB = {psat:.9f}; equilibrium states: {saturation}')


if __name__=='__main__':main()
