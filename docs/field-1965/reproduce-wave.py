"""Recompute Field (1965), Figure 3: the complex wave pair, not the fastest root.

Run: python reproduce-wave.py --output output
Keep reproduce.py in the same folder; it supplies the cubic and fluid matrix.
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

from reproduce import linear_eigenvalues, polynomial


def wave_root(q, alpha, beta, gamma):
    """Select the root with positive imaginary part (exp(+i omega t)).

    Fail explicitly if the complex pair ceases to be unique. No sorting by
    growth rate or continuity heuristic can silently substitute condensation.
    At q=0 use the wave-pair limit x=0; the third root is -alpha.
    """
    if q < 0 or beta < 0 or gamma <= 1:
        raise ValueError("Requires q>=0, beta>=0, gamma>1")
    if q == 0:
        return 0j
    s = alpha + beta*q*q
    roots = np.roots([1.0, s, q*q, q*q*(s-1)/gamma])
    pair = [complex(r) for r in roots if r.imag > 1e-12*max(1,abs(r))]
    if len(pair) != 1:
        raise ValueError(f"No unique oscillatory wave pair at q={q}, beta={beta}")
    return pair[0]


def real_part_equation(xi, q, alpha, beta, gamma):
    """Equation (40), regrouped using s=alpha+beta*q^2 and xi=Re(x)."""
    s = alpha + beta*q*q
    return xi**3 + s*xi**2 + (s*s+q*q)*xi/4 + q*q*((1-1/gamma)*s+1/gamma)/8


def validate(config):
    a,g = config["alpha"],config["gamma"]
    if not (a == -2 and np.isclose(g,5/3) and config["betas"] == [0,.01,.1,1]):
        raise ValueError("This validated exhibit is restricted to the Figure 3 parameters")
    checks = []
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    def condition(name,passed,**measurements):
        checks.append(dict(name=name,passed=bool(passed),**measurements))
    eigen_errors,cubic_errors,real_errors,conjugate_errors=[],[],[],[]
    for b in config["betas"]:
        for q in np.geomspace(1e-5,config["q_max"],240):
            root=wave_root(q,a,b,g)
            ev=linear_eigenvalues(q,a,b,g)
            eigen_errors.append(min(abs(root-r) for r in ev)/max(1,abs(root)))
            conjugate_errors.append(min(abs(root.conjugate()-r) for r in ev)/max(1,abs(root)))
            s=a+b*q*q
            scale=1+abs(root**3)+abs(s*root*root)+abs(q*q*root)+abs(q*q*(s-1)/g)
            cubic_errors.append(abs(polynomial(root,q,a,b,g))/scale)
            xi=root.real
            scale40=1+abs(xi**3)+abs(s*xi*xi)+abs((s*s+q*q)*xi/4)+abs(q*q*((1-1/g)*s+1/g)/8)
            real_errors.append(abs(real_part_equation(xi,q,a,b,g))/scale40)
    check("maximum scaled wave-root error against fluid-matrix eigenvalues",max(eigen_errors),1e-10)
    check("conjugate wave root also matches fluid-matrix eigenvalues",max(conjugate_errors),1e-10)
    check("maximum scaled equation (39) complex residual",max(cubic_errors),1e-12)
    check("maximum scaled equation (40) real-part residual",max(real_errors),1e-12)
    q=1e-3
    expected_growth=-(1+(g-1)*a)/(2*g*a*a)
    for b in config["betas"]:
        root=wave_root(q,a,b,g)
        check(f"long-wave Re(x)/q^2 limit, beta={b:g}",abs(root.real/q**2-expected_growth),1e-6)
    check("long-wave Im(x)/q limit",abs(wave_root(q,a,0,g).imag/q-np.sqrt((a-1)/(g*a))),1e-6)
    limit=-(1+(g-1)*a)/(2*g)
    check("zero-conduction short-wave Re(x) limit",abs(wave_root(1e5,a,0,g).real-limit),1e-8)
    check("zero-conduction short-wave Im(x)/q limit",abs(wave_root(1e5,a,0,g).imag/1e5-1),1e-8)
    for b in config["betas"]:
        if not b:continue
        qc=float(np.sqrt((-a-1/(g-1))/b))
        at=wave_root(qc,a,b,g)
        check(f"neutral Re(x) at wave cutoff, beta={b:g}",abs(at.real),1e-12)
        check(f"Im(x)=q at wave cutoff, beta={b:g}",abs(at.imag-qc),1e-12)
        before=wave_root(qc*(1-1e-4),a,b,g).real
        after=wave_root(qc*(1+1e-4),a,b,g).real
        condition(f"wave changes from growth to damping, beta={b:g}",before>0 and after<0,
                  below=float(before),above=float(after))
        fastest=float(max(linear_eigenvalues(qc*(1+1e-4),a,b,g).real))
        condition(f"condensation still grows above wave cutoff, beta={b:g}",fastest>0,
                  fastest_root_real=fastest)
    root=wave_root(.1,a,0,g)
    fastest=float(max(linear_eigenvalues(.1,a,0,g).real))
    condition("wave selector does not select fastest condensation root",root.real<.01 and fastest>1,
              wave_growth=float(root.real),fastest_root_real=fastest)
    if not all(c["passed"] for c in checks):
        raise RuntimeError(json.dumps(checks,indent=2))
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path(__file__).with_name("config-wave.json"))
    parser.add_argument("--output",type=Path,default=Path("output"))
    args=parser.parse_args()
    config=json.loads(args.config.read_text())
    checks=validate(config)
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    a,g=config["alpha"],config["gamma"]
    fig,ax=plt.subplots(figsize=(9,5.8),layout="constrained")
    colors=["#263735","#22645c","#697570","#8b7050"]
    styles=["-","--","-.",":"]
    rows,peaks=[],[]
    for b,color,style in zip(config["betas"],colors,styles):
        qc=float(np.sqrt((-a-1/(g-1))/b)) if b else None
        q=np.linspace(0,config["q_max"],config["samples"])
        if qc is not None and qc<=config["q_max"]:q=np.unique(np.append(q,qc))
        roots=np.array([wave_root(float(qi),a,b,g) for qi in q])
        plotted=roots.real.copy()
        if qc is not None:
            plotted[q>qc]=np.nan
            plotted[np.isclose(q,qc,rtol=0,atol=1e-14)]=0
        ax.plot(q,plotted,label=rf"$\beta={b:g}$",color=color,ls=style,lw=2)
        for qi,root in zip(q,roots):
            neutral=qi==0 or (qc is not None and abs(qi-qc)<1e-14)
            domain="neutral" if neutral else "growing_wave" if qc is None or qi<qc else "damped_wave"
            rows.append([b,float(qi),float(root.real),float(root.imag),domain])
        idx=int(np.nanargmax(plotted))
        peaks.append(dict(beta=b,wave_cutoff_q=qc,sampled_peak_q=float(q[idx]),
                          sampled_peak_real_x=float(plotted[idx]),peak_is_interior=bool(idx<len(q)-1)))
    limit=-(1+(g-1)*a)/(2*g)
    ax.axhline(limit,color="#888888",lw=1,ls=(0,(3,3)))
    ax.text(6.95,.101,"Zero-conduction short-wave limit: 0.1",ha="right",va="bottom",fontsize=9,color="#555555")
    ax.set(xlim=(0,7),ylim=(0,.108),xlabel=r"Normalized wavenumber $k/k_\rho$ (dimensionless)",
           ylabel=r"Wave growth rate $\mathrm{Re}(n)/(k_\rho c)$ (dimensionless)",
           title="Field (1965), Figure 3 — wave mode")
    ax.text(.02,.96,r"$\alpha=-2,\quad \gamma=5/3$",transform=ax.transAxes,va="top",fontsize=12)
    ax.legend(loc="center right",frameon=False,title="Thermal conduction")
    ax.spines[["top","right"]].set_visible(False);ax.grid(alpha=.15)
    fig.supxlabel("Recomputed from Eq. (39); growing wave branches and neutral endpoints are shown.\nSource: Field, ApJ 142, 531–567 (1965), pp. 537–540; the condensation root is excluded.",fontsize=9)
    fig.savefig(out/"figure-3.png",dpi=180)
    fig.savefig(out/"figure-3.svg",metadata={"Date":None});plt.close(fig)
    with (out/"figure-3.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["beta","k_over_k_rho","wave_real_n_over_k_rho_c","wave_imag_n_over_k_rho_c","domain"]);w.writerows(rows)
    files={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
           for name in ["reproduce-wave.py","reproduce.py"]}
    report=dict(source="Field (1965), Figure 3; equations (15), (27), (30), (38)-(40)",
                config=config,checks=checks,curves=peaks,
                limitations=["No digitized pointwise comparison with the printed curves.",
                             "Modern line styles and typography; not a facsimile.",
                             "Plot omits damped wave branches; CSV retains their signed real and positive imaginary parts.",
                             "Wave stabilization is not full-system stabilization: condensation can remain unstable.",
                             "Peaks are sampled maxima; beta=0 has no finite-wavenumber maximum."],
                environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
                code_sha256=files,config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/"validation-wave.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(checks_passed=len(checks),output=str(out),curves=peaks),indent=2))


if __name__=="__main__":main()
