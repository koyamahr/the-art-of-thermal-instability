"""FGH/GHF 1969 model ingredients, not a reproduction of equilibrium figures.

Run: python reproduce.py --output output
Inputs: config.json (visually transcribed GHF Table 1) and source-audit.json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save_figure(fig, out, name):
    fig.savefig(out/(name+'.png'),dpi=180)
    svg=out/(name+'.svg')
    fig.savefig(svg,metadata={'Date':None})
    # Keep generated text friendly to Git's whitespace checks.
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')


def heating(x):
    """GHF (1969) equations (23), (25), and (26), preserving rounded prefactors."""
    x=np.asarray(x,dtype=float)
    if np.any(x<0):raise ValueError("Electron ratio cannot be negative")
    ratio=(1+130*x)/(1+57*x)
    return 9*ratio, .13*ratio, 1.4e-11*ratio


def log_coefficient(temperature,transition,collider):
    """Natural logarithm of an energy-loss coefficient, equation (16).

    Coefficient units: erg cm^3 s^-1. R carries an additional K^-m.
    Keep the log representation to preserve extremely small low-T values.
    """
    t=np.asarray(temperature,dtype=float)
    if np.any(t<=0):raise ValueError("Temperature must be positive")
    fit=transition[collider]
    if fit is None:raise ValueError("No coefficient tabulated for this collider")
    return np.log(fit["R"])+fit["m"]*np.log(t)-transition["T_star_K"]/t


def validate(config):
    checks=[]
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,passed,**details):checks.append(dict(name=name,passed=bool(passed),**details))
    # Published approximate anchors in the paragraph following equation (26).
    for x,quoted in [(1e-4,9),(.01,13),(1,20)]:
        check(f"equation (26) agrees with rounded {quoted} eV anchor at x={x:g}",abs(float(heating(x)[0])-quoted),.5)
    check("zero-electron limit of equation (26)",abs(float(heating(0)[0])-9),1e-12)
    check("large-electron limit of equation (26)",abs(float(heating(1e10)[0])-9*130/57),1e-9)
    x=np.geomspace(1e-4,1,501);energy,f,normalization=heating(x)
    condition("heating increases while remaining within the 20 eV input energy to printed rounding",np.all(np.diff(energy)>0) and np.max(energy)<20.5,
              maximum_printed_fit_eV=float(np.max(energy)),note="The rounded fit slightly exceeds 20 eV at x=1; it is not extrapolated as a physical bound.")
    condition("efficiency stays between zero and one",np.all((f>0)&(f<1)))
    # Independent centered finite difference checks the temperature exponent and sign.
    derivative_errors=[];pairs=0
    for row in config["transitions"]:
        for collider in ["electron","hydrogen"]:
            if row[collider] is None:continue
            pairs+=1
            for t in [100,1000,10000]:
                h=1e-4
                measured=(log_coefficient(t*np.exp(h),row,collider)-log_coefficient(t*np.exp(-h),row,collider))/(2*h)
                expected=row[collider]["m"]+row["T_star_K"]/t
                derivative_errors.append(abs(measured-expected)/max(1,abs(expected)))
    check("logarithmic cooling slopes match m+T_star/T",max(derivative_errors),1e-7)
    condition("Table 1 contains five electron and six H-atom fits",pairs==11 and
              sum(r['electron'] is not None for r in config['transitions'])==5 and
              sum(r['hydrogen'] is not None for r in config['transitions'])==6,
              electron_fits=sum(r['electron'] is not None for r in config['transitions']),
              hydrogen_fits=sum(r['hydrogen'] is not None for r in config['transitions']))
    rejected=0
    for row in config['transitions']:
        for collider in ['electron','hydrogen']:
            if row[collider] is None:
                try:log_coefficient(100,row,collider)
                except ValueError:rejected+=1
    condition("blank coefficients are unavailable rather than zero",rejected==3,unavailable_pairs=rejected)
    # Physically independent limiting cross-over: same activation temperature
    # cancels when comparing collider contributions for a single transition.
    row=config['transitions'][0];t=100
    critical_ratio=np.exp(log_coefficient(t,row,'hydrogen')-log_coefficient(t,row,'electron'))
    reference=(1.7e-23/7.9e-20)*np.sqrt(t)
    check("C II equal-contribution n_e/n_H ratio at 100 K",abs(float(critical_ratio/reference)-1),1e-12)
    if not all(c['passed'] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('config.json'))
    parser.add_argument('--audit',type=Path,default=Path(__file__).with_name('source-audit.json'))
    parser.add_argument('--output',type=Path,default=Path('output'))
    args=parser.parse_args();config=json.loads(args.config.read_text());audit=json.loads(args.audit.read_text())
    checks=validate(config);out=args.output;out.mkdir(parents=True,exist_ok=True)
    x=np.geomspace(config['electron_ratio_min'],config['electron_ratio_max'],config['samples'])
    energy,f,norm=heating(x)
    fig,axes=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
    axes[0].semilogx(x,energy,color='#22645c',lw=2)
    anchors=np.array([1e-4,.01,1]);axes[0].scatter(anchors,heating(anchors)[0],color='#22645c',s=25)
    axes[0].set(ylabel='Indirect heat per primary electron (eV)',ylim=(8,21),title='Energy transferred to the gas · Eq. (26)')
    axes[1].semilogx(x,f,color='#8b7050',lw=2)
    axes[1].set(ylabel='Heating efficiency f (dimensionless)',ylim=(.1,.32),title='Fraction of ionization + excitation loss · Eq. (25)')
    for ax in axes:
        ax.set_xlabel(r'Electron ratio $x=n_e/n(\mathrm{H\ I})$ (dimensionless)')
        ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Cosmic-ray heating · ingredients of the FGH/GHF 1969 model',fontsize=14)
    fig.supxlabel('Derived diagnostics from Goldsmith, Habing & Field, ApJ 158, 173 (1969), p. 181.\nThese are evaluations of printed fits, not the equilibrium figures of the FGH Letter.',fontsize=9)
    save_figure(fig,out,'heating')
    plt.close(fig)
    with (out/'heating.csv').open('w',newline='') as file:
        w=csv.writer(file,lineterminator='\n');w.writerow(['electron_ratio_ne_over_nHI','indirect_heat_eV_eq26','efficiency_eq25','Gamma_ion_over_zeta_n_target_erg_eq23']);w.writerows(zip(x,energy,f,norm))
    temperatures=np.geomspace(config['temperature_min_K'],config['temperature_max_K'],config['samples'])
    fig,axes=plt.subplots(1,2,figsize=(10,5.7),layout='constrained')
    colors=['#22645c','#8b7050','#596662','#596662','#776b80','#776b80','#776b80']
    styles=['-','-','-','--','-','--',':'];rows=[]
    for ax,collider,title in zip(axes,['electron','hydrogen'],['Electron excitation · 5 tabulated fits','H-atom excitation · 6 tabulated fits']):
        for row,color,style in zip(config['transitions'],colors,styles):
            if row[collider] is None:continue
            log_r=log_coefficient(temperatures,row,collider)
            ax.loglog(temperatures,np.exp(log_r),color=color,ls=style,lw=1.6,label=row['label'])
            rows.extend([row['id'],collider,float(t),float(lr/np.log(10))] for t,lr in zip(temperatures,log_r))
        ax.set(xlim=(10,20000),ylim=(1e-29,1e-18),xlabel='Temperature T (K)',ylabel=r'Energy-loss coefficient $r_k$ (erg cm$^3$ s$^{-1}$)',title=title)
        ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8,loc='upper left');ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Metal cooling coefficients · GHF (1969), Table 1 and Eq. (16)',fontsize=14)
    fig.supxlabel('Fits evaluated over 10–20,000 K; display clips coefficients below 10⁻²⁹. CSV retains log₁₀ coefficients.\nMissing Table 1 entries are omitted. These coefficients alone do not determine equilibrium or cooling dominance.',fontsize=9)
    save_figure(fig,out,'cooling-coefficients')
    plt.close(fig)
    with (out/'cooling-coefficients.csv').open('w',newline='') as file:
        w=csv.writer(file,lineterminator='\n');w.writerow(['transition','collider','temperature_K','log10_r_erg_cm3_s']);w.writerows(rows)
    hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['reproduce.py','source-audit.json']}
    report=dict(scope=audit['scope'],source=config['source'],checks=checks,
                environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
                code_sha256=hashes,config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest(),
                limitations=['Derived diagnostics, not a reproduction of FGH equilibrium Figures 1 or 2.',
                             'This ingredient script does not solve equilibrium; cold_branch.py implements a separate reduced model.',
                             'Printed rounded heating fits are preserved; their normalizations are not forced to agree exactly.',
                             'Temperature evaluation range is a diagnostic choice, not an independently established validity range of every historical fit.',
                             'Cooling-curve weights require actual populations, densities, and abundances.'])
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(checks_passed=len(checks),tabulated_collider_fits=11,scope=report['scope'],output=str(out)),indent=2))


if __name__=='__main__':main()
