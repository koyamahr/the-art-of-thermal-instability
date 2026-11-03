"""Recompute Field (1965), Figure 2, from equations (15), (27), (38), (39).

Run: python reproduce.py --output output
No source PDF or external data is required to run the calculation.
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


def polynomial(x, q, alpha, beta, gamma):
    """Equation (39), regrouped as a cubic in x; q = k/k_rho."""
    u = q * q
    return ((x + alpha + beta * u) * x + u) * x + u * (alpha - 1 + beta * u) / gamma


def growth(q, alpha, beta, gamma):
    """Unique nonnegative condensation root in this Figure 2 parameter regime.

    Return NaN in the stable domain rather than silently plotting an acoustic
    branch. At q=0 return the limiting value 0 (the cubic has a double zero).
    The function is deliberately restricted to 0 < alpha < 1, beta >= 0.
    """
    if not (0 < alpha < 1 and beta >= 0 and gamma > 1 and q >= 0):
        raise ValueError("Requires 0 < alpha < 1, beta >= 0, gamma > 1, q >= 0")
    if q == 0:
        return 0.0
    if beta > 0:
        cutoff = np.sqrt((1 - alpha) / beta)
        if np.isclose(q, cutoff, rtol=0, atol=1e-14):
            return 0.0
        if q > cutoff:
            return float("nan")
    # P(0)<0 and P((1-alpha)/gamma)>0 inside the unstable band.
    # P'(x)>0 for x>=0: bisection cannot switch to another mode.
    low, high = 0.0, (1 - alpha) / gamma
    for _ in range(100):
        mid = (low + high) / 2
        if polynomial(mid, q, alpha, beta, gamma) > 0:
            high = mid
        else:
            low = mid
    return float((low + high) / 2)


def linear_eigenvalues(q, alpha, beta, gamma):
    """Independent linear-system check using equations (11)-(14).

    State: (delta rho/rho_0, i delta v/c, delta p/p_0).
    Dimensionless time: tau=k_rho*c*t. Equation of state gives
    delta T/T_0 = delta p/p_0 - delta rho/rho_0.
    """
    s = alpha + beta * q * q
    matrix = np.array([[0, -q, 0], [0, 0, q / gamma],
                       [s - 1, -gamma * q, -s]], dtype=float)
    return np.linalg.eigvals(matrix)


def validate(config):
    a, g = config["alpha"], config["gamma"]
    checks = []
    def check(name, measured, tolerance):
        checks.append(dict(name=name, error=float(measured), tolerance=tolerance,
                           passed=bool(measured < tolerance)))
    errors, residuals = [], []
    for b in config["betas"]:
        upper = min(config["q_max"], np.sqrt((1-a)/b) * (1-1e-7)) if b else config["q_max"]
        for q in np.geomspace(1e-5, upper, 160):
            x = growth(q, a, b, g)
            eigenvalues = linear_eigenvalues(q, a, b, g)
            errors.append(abs(x - max(eigenvalues.real)))
            u = q*q
            scale = 1 + abs(x**3) + abs((a+b*u)*x*x) + abs(u*x) + abs(u*(a-1+b*u)/g)
            residuals.append(abs(polynomial(x,q,a,b,g))/scale)
    check("maximum absolute error against linear-system eigenvalues", max(errors), 1e-10)
    check("maximum scaled cubic residual", max(residuals), 1e-12)
    q = 1e-7
    check("long-wave x/q limit, beta=0", abs(growth(q,a,0,g)/q - np.sqrt((1-a)/(g*a))), 1e-5)
    check("short-wave x limit, beta=0", abs(growth(1e5,a,0,g) - (1-a)/g), 1e-9)
    for b in config["betas"]:
        if not b:
            continue
        qc = float(np.sqrt((1-a)/b))
        check(f"neutral root at cutoff, beta={b:g}", abs(polynomial(0,qc,a,b,g)), 1e-12)
        before = max(linear_eigenvalues(qc*(1-1e-4),a,b,g).real)
        after = max(linear_eigenvalues(qc*(1+1e-4),a,b,g).real)
        checks.append(dict(name=f"growth changes to damping across cutoff, beta={b:g}",
                           below=float(before), above=float(after), passed=bool(before>0 and after<0)))
    assert all(c["passed"] for c in checks), json.dumps(checks, indent=2)
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.json"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    checks = validate(config)
    a,g = config["alpha"],config["gamma"]
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9,5.8), layout="constrained")
    colors = ["#263735", "#22645c", "#697570", "#8b7050"]
    styles = ["-", "--", "-.", ":"]
    rows, peaks = [], []
    for b, color, style in zip(config["betas"], colors, styles):
        qc = float(np.sqrt((1-a)/b)) if b else None
        q = np.linspace(0, config["q_max"], config["samples"])
        if qc is not None and qc <= config["q_max"]:
            q = np.unique(np.append(q,qc))
        x = np.array([growth(float(v),a,b,g) for v in q])
        ax.plot(q,x,label=rf"$\beta={b:g}$",color=color,ls=style,lw=2)
        for qi,xi in zip(q,x):
            rows.append([b,float(qi),None if np.isnan(xi) else float(xi),
                         "stable_omitted" if np.isnan(xi) else "neutral" if xi==0 else "unstable"])
        idx=int(np.nanargmax(x))
        peaks.append(dict(beta=b,cutoff_q=qc,sampled_peak_q=float(q[idx]),
                          sampled_peak_x=float(x[idx]),peak_is_interior=bool(idx<len(q)-1)))
    ax.axhline((1-a)/g,color="#888888",lw=1,ls=(0,(3,3)))
    ax.text(6.95,.303,"Zero-conduction short-wave limit: 0.3",ha="right",va="bottom",fontsize=9,color="#555555")
    ax.set(xlim=(0,7),ylim=(0,.325),xlabel=r"Normalized wavenumber $k/k_\rho$ (dimensionless)",
           ylabel=r"Normalized growth rate $n/(k_\rho c)$ (dimensionless)",
           title="Field (1965), Figure 2 — condensation mode")
    ax.text(.02,.96,r"$\alpha=1/2,\quad \gamma=5/3$",transform=ax.transAxes,va="top",fontsize=12)
    ax.legend(loc="center right",frameon=False,title="Thermal conduction")
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(alpha=.15)
    fig.supxlabel("Recomputed from Eq. (39); only the unstable branch and its neutral endpoints are shown.\nSource: Field, ApJ 142, 531–567 (1965), pp. 534, 537, 539–540.",fontsize=9)
    fig.savefig(out/"figure-2.png",dpi=180)
    fig.savefig(out/"figure-2.svg",metadata={"Date":None})
    plt.close(fig)
    with (out/"figure-2.csv").open("w",newline="") as f:
        writer=csv.writer(f);writer.writerow(["beta","k_over_k_rho","n_over_k_rho_c","domain"]);writer.writerows(rows)
    report=dict(source="Field (1965), Figure 2; equations (11)-(16), (27), (38)-(39)",
                config=config,checks=checks,curves=peaks,
                limitations=["Analytical and numerical checks; no digitized pointwise comparison with the printed curves.",
                             "Modern typography and line styles; not a facsimile of the original figure.",
                             "Only the unstable condensation branch is drawn. Stable-domain CSV growth entries are blank.",
                             "The beta=0 sampled maximum is a finite-domain endpoint, not a physical finite-wavenumber maximum."],
                environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                config_sha256=hashlib.sha256(args.config.read_bytes()).hexdigest())
    (out/"validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(checks_passed=len(checks),output=str(out),curves=peaks),indent=2))


if __name__ == "__main__":
    main()
