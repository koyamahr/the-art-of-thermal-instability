"""Recompute Field (1965), Figure 1, from equations (17)-(28).

Run: python reproduce-map.py --output output
No original PDF or helper module is required at runtime.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def roots(t, r, gamma):
    """Equation (18): y=n/(k*c), t=sigma_T_prime, r=sigma_rho."""
    return np.roots([1, t, 1, (t-r)/gamma])


def discriminant(t, r, gamma):
    """Field's Cardano D, equation (19), NOT the conventional discriminant.

    D<0: three real distinct roots. D>0: one real and a conjugate pair.
    The conventional cubic discriminant is -108*D.
    """
    c=(t-r)/gamma
    return (c-t/3+2*t**3/27)**2/4 + (1-t*t/3)**3/27


def positive_count(t, r, gamma):
    """Count roots with Re(y)>0 away from neutral lines, including multiplicity."""
    return np.asarray(r>t,dtype=int)+2*np.asarray(r<-(gamma-1)*t,dtype=int)


def thorn(t, gamma):
    """Equation (20): lower and upper r boundaries for |t|>=sqrt(3)."""
    t=np.asarray(t)
    center=(1-gamma/3)*t+2*gamma*t**3/27
    delta=2*gamma*np.maximum(t*t-3,0)**1.5/27
    return center-delta,center+delta


def matched_error(a,b):
    """Root-set distance independent of eigensolver ordering."""
    return min(max(abs(x-y) for x,y in zip(a,p)) for p in itertools.permutations(b))


def validate(config):
    g=config["gamma"];checks=[]
    if not (np.isclose(g,5/3) and config["alphas"]==[.5,-2] and config["betas"]==[.01,.1,1]):
        raise ValueError("This exhibit requires the original Figure 1 parameters")
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,passed=bool(error<tolerance)))
    count_mismatches=0;topology_mismatches=0;eigen_errors=[];residuals=[];tested=0
    for t in np.linspace(-7.93,7.93,53):
        for r in np.linspace(-7.89,7.89,51):
            if min(abs(t-r),abs((g-1)*t+r),abs(discriminant(t,r,g)))<1e-6:continue
            rr=roots(t,r,g)
            count_mismatches+=int(np.count_nonzero(rr.real>1e-8)!=positive_count(t,r,g))
            all_real=np.max(abs(rr.imag))<1e-8
            topology_mismatches+=int(all_real!=(discriminant(t,r,g)<0))
            matrix=np.array([[0,-1,0],[0,0,1/g],[t-r,-g,-t]])
            eigen_errors.append(matched_error(rr,np.linalg.eigvals(matrix)))
            for y in rr:
                residuals.append(abs(y**3+t*y*y+y+(t-r)/g)/(1+abs(y**3)+abs(t*y*y)+abs(y)+abs((t-r)/g)))
            tested+=1
    check("positive-root count mismatches on off-boundary grid",count_mismatches,.5)
    check("Cardano D topology mismatches on off-boundary grid",topology_mismatches,.5)
    check("maximum root-set error against fluid-matrix eigenvalues",max(eigen_errors),1e-10)
    check("maximum scaled equation (18) residual",max(residuals),1e-12)
    checks[0]["points_tested"]=tested;checks[1]["points_tested"]=tested
    boundary_errors=[]
    for sign in [-1,1]:
        t=sign*np.linspace(np.sqrt(3),8,300)
        for r in thorn(t,g):
            for ti,ri in zip(t,r):
                c=(ti-ri)/g
                scale=1+abs((c-ti/3+2*ti**3/27)**2/4)+abs((1-ti*ti/3)**3/27)
                boundary_errors.append(abs(discriminant(ti,ri,g))/scale)
    check("equation (20) boundaries satisfy D=0",max(boundary_errors),1e-12)
    for sign in [-1,1]:
        t=sign*np.sqrt(3);r=(1-g/9)*t;y=-t/3
        check(f"triple-root polynomial and derivatives at cusp sign={sign}",
              max(abs(y**3+t*y*y+y+(t-r)/g),abs(3*y*y+2*t*y+1),abs(6*y+2*t)),1e-12)
    neutral_errors=[]
    for t in np.linspace(-8,8,61):
        neutral_errors.append(abs((t-t)/g))
        # At r=-(gamma-1)t the polynomial factors (y+t)(y^2+1).
        r=-(g-1)*t
        for y in [1j,-1j,-t]:neutral_errors.append(abs(y**3+t*y*y+y+(t-r)/g))
    check("neutral-line zero and imaginary roots",max(neutral_errors),1e-12)
    path_errors=[]
    for a in config["alphas"]:
        for b in config["betas"]:
            for q in np.geomspace(.15,7,50):
                r=1/q;t=a/q+b*q
                xroots=np.roots([1,a+b*q*q,q*q,q*q*(a-1+b*q*q)/g])
                path_errors.append(matched_error(q*roots(t,r,g),xroots))
    check("path mapping y=x/q agrees with Figures 2 and 3 cubic",max(path_errors),1e-10)
    for b in config["betas"]:
        a=.5;q=np.sqrt((1-a)/b)
        check(f"Figure 2 path crosses isobaric neutral line, beta={b:g}",abs(a/q+b*q-1/q),1e-12)
        a=-2;q=np.sqrt((-a-1/(g-1))/b)
        check(f"Figure 3 path crosses isentropic neutral line, beta={b:g}",abs((g-1)*(a/q+b*q)+1/q),1e-12)
    if not all(c["passed"] for c in checks):raise RuntimeError(json.dumps(checks,indent=2))
    return checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path(__file__).with_name("config-map.json"))
    parser.add_argument("--output",type=Path,default=Path("output"))
    args=parser.parse_args();config=json.loads(args.config.read_text());checks=validate(config)
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    g=config["gamma"];limit=config["extent"]
    axis=np.linspace(-limit,limit,config["grid_samples"])
    tt,rr=np.meshgrid(axis,axis);counts=positive_count(tt,rr,g);dd=discriminant(tt,rr,g)
    fig,ax=plt.subplots(figsize=(9,11))
    fig.subplots_adjust(left=.12,right=.98,bottom=.26,top=.96)
    ax.contourf(tt,rr,counts,levels=[-.5,.5,1.5,2.5,3.5],colors=["#f5f4ee","#e3ebe7","#ece9df","#d3ddd6"])
    for sign in [-1,1]:
        t=sign*np.linspace(np.sqrt(3),limit,1600)
        lo,hi=thorn(t,g)
        ax.fill_between(t,lo,hi,facecolor="none",hatch="///",edgecolor="#83918a",linewidth=0,zorder=2)
        ax.plot(t,lo,color="#55645d",lw=1.2,zorder=3)
        ax.plot(t,hi,color="#55645d",lw=1.2,zorder=3)
    ax.plot(axis,axis,color="#263735",lw=1.6,zorder=4)
    ax.plot(axis,-(g-1)*axis,color="#263735",ls="--",lw=1.6,zorder=4)
    for a,color in zip(config["alphas"],["#22645c","#8b7050"]):
        rho=np.geomspace(.001,limit,config["path_samples"])
        ax.plot(a*rho,rho,color=color,ls=(0,(1,3)),lw=1,alpha=.8)
        for b,style in zip(config["betas"],["-","--","-."]):
            t=a*rho+b/rho
            ax.plot(t,rho,color=color,ls=style,lw=1.5,zorder=5)
    bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=4)
    for x,y,label in [(5.3,.5,"0 growing roots"),(.0,6.1,"1 growing root"),(-5.3,-.5,"3 growing roots"),(.0,-6.1,"2 growing roots")]:
        ax.text(x,y,label,ha="center",va="center",fontsize=12,bbox=bbox,zorder=6)
    ax.text(5.2,6.9,"Hatching: 3 real roots",ha="center",fontsize=10,bbox=bbox,zorder=6)
    ax.annotate("Triple root",xy=(np.sqrt(3),(1-g/9)*np.sqrt(3)),xytext=(3.6,2.1),
                fontsize=10,arrowprops=dict(arrowstyle="->",color="#55645d"),bbox=bbox,zorder=6)
    ax.axhline(0,color="#a5aca5",lw=.7);ax.axvline(0,color="#a5aca5",lw=.7)
    ax.set(xlim=(-limit,limit),ylim=(-limit,limit),aspect="equal",
           xlabel=r"$\sigma_T^\prime=(k_T+k^2/k_K)/k$ (dimensionless)",
           ylabel=r"$\sigma_\rho=k_\rho/k$ (dimensionless)",title="Field (1965), Figure 1 — domains of instability")
    ax.text(.02,.97,r"$\gamma=5/3$",transform=ax.transAxes,va="top",fontsize=12,bbox=bbox)
    handles=[Line2D([],[],color="#263735",label=r"Isobaric neutral: $\sigma_\rho=\sigma_T^\prime$"),
             Line2D([],[],color="#263735",ls="--",label=r"Isentropic neutral: $\sigma_\rho=-\frac{2}{3}\sigma_T^\prime$"),
             Line2D([],[],color="#22645c",label=r"Paths: $\alpha=1/2$ (Figure 2)"),
             Line2D([],[],color="#8b7050",label=r"Paths: $\alpha=-2$ (Figure 3)"),
             Line2D([],[],color="#666666",ls="-",label=r"$\beta=0.01$"),
             Line2D([],[],color="#666666",ls="--",label=r"$\beta=0.1$"),
             Line2D([],[],color="#666666",ls="-.",label=r"$\beta=1$"),
             Line2D([],[],color="#666666",ls=(0,(1,3)),label=r"$\beta=0$ reference rays")]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.55,.055),ncols=2,frameon=False,fontsize=9)
    fig.text(.5,.015,"Recomputed from Eqs. (18)–(20), (25), (28); growing means Re(y)>0, counting multiplicity.\nSource: Field, ApJ 142, 531–567 (1965), pp. 535–537. Neutral boundaries are treated separately.",fontsize=8.5,ha="center",va="bottom")
    fig.savefig(out/"figure-1.png",dpi=180)
    fig.savefig(out/"figure-1.svg",metadata={"Date":None});plt.close(fig)
    with (out/"figure-1-grid.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["sigma_T_prime","sigma_rho","cardano_D","growing_root_count","root_topology"])
        for t,r,d,count in zip(tt.flat,rr.flat,dd.flat,counts.flat):
            neutral=min(abs(t-r),abs((g-1)*t+r))<1e-10
            topology="repeated_root_boundary" if abs(d)<1e-10 else "three_real" if d<0 else "one_real_two_complex"
            w.writerow([float(t),float(r),float(d),"" if neutral else int(count),topology])
    with (out/"figure-1-paths.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["alpha","beta","k_over_k_rho","sigma_T_prime","sigma_rho"])
        for a in config["alphas"]:
            for b in config["betas"]:
                for rho in np.geomspace(.001,limit,config["path_samples"]):
                    t=a*rho+b/rho
                    if abs(t)<=limit:w.writerow([a,b,float(1/rho),float(t),float(rho)])
    report=dict(source="Field (1965), Figure 1; equations (17)-(28)",config=config,checks=checks,
                limitations=["No digitized pointwise comparison with the printed figure.",
                             "Root-count background, typography, and legend are modernized.",
                             "The original circular annotation guide is omitted: it is not a stability boundary.",
                             "Neutral-line grid counts are blank; topology labels use an absolute 1e-10 D tolerance.",
                             "Grid root counts include multiplicity and count Re(y)>0, not merely real roots."],
                environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/"validation-map.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(checks_passed=len(checks),grid_points=tt.size,output=str(out)),indent=2))


if __name__=="__main__":main()
