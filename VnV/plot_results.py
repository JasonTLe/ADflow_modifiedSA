"""Plot ADflow results of the validation cases against the experimental data.

Usage: MPLBACKEND=Agg python plot_results.py [plate|fin|flare|all]
Reads results/adflow/<case>/medium/<model>/lines.json (fin, flare) or
cells.npz + bl.json (plate) and data/experiments_nasa.json (tabulated NASA
data, see data/README.md) and writes results/figures/<case>_vs_experiment.png."""

import os
import sys
import json
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
from common import FIGURES, run_dir  # noqa: E402
EXP = json.load(open(os.path.join(HERE, "data", "experiments_nasa.json")))
STYLE = {
    "SA": dict(color="k", ls="-"),
    "SA-Edwards": dict(color="tab:blue", ls="-"),
    "SA-comp": dict(color="tab:purple", ls="--"),
    "SA-Catris": dict(color="tab:red", ls="-"),
    "SA-Catris-cons": dict(color="tab:orange", ls="--"),
    "SA-Catris1999": dict(color="tab:brown", ls=":"),
    "SA-Seror": dict(color="tab:green", ls="-."),
    "SA-cons": dict(color="tab:gray", ls="-."),
    "SA-Edwards-Catris": dict(color="tab:cyan", ls="--"),
}
LABEL = {"SA": "SA (baseline)", "SA-Edwards": "SA-Edwards (Edwards & Chandra)", "SA-comp": "SA-comp (Spalart compressibility)",
         "SA-Catris": "SA-Catris (Catris-Aupoix)", "SA-Catris-cons": "SA-Catris, conservative form",
         "SA-Catris1999": "SA-Catris (1999 form)", "SA-Seror": "SA-Seror", "SA-cons": "SA, conservative diffusion",
         "SA-Edwards-Catris": "SA-Edwards + Catris"}
EXPSTYLE = dict(marker="o", mfc="none", color="k", ls="none", ms=5)
EXPSTYLE2 = dict(marker="s", mfc="none", color="tab:gray", ls="none", ms=4)


def load(case, fileName="lines.json"):
    out = {}
    for f in sorted(glob.glob(os.path.join(run_dir("adflow", case, "medium", "*"), fileName))):
        d = json.load(open(f))
        name = os.path.basename(os.path.dirname(f))
        if name in STYLE:
            out[name] = d
    return {m: out[m] for m in STYLE if m in out}


def savefig(fig, name):
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, name), dpi=130)


def plot_plate():
    """Incoming boundary layer of the Kussoy & Horstman experiment (survey at
    x = -11 cm) vs the 2-D plate solutions: profiles of u/U, T/Tinf and M,
    and the integral/wall quantities."""
    e = EXP["kussoy_fin15"]
    prof = e["upstream_BL_profile"]
    bl = load("plate", "bl.json")
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.8))
    for m in bl:
        f = os.path.join(run_dir("adflow", "plate", "medium", m), "cells.npz")
        if not os.path.exists(f):
            continue
        d = np.load(f)
        x, w = d["x"], d["w"]
        xu = np.unique(np.round(x[:, 0], 8))
        xa = xu[np.argmin(np.abs(xu + 0.11))]
        sel = (np.abs(x[:, 0] - xa) < 1e-7) & (x[:, 2] < x[:, 2].min() + 1e-9)
        o = np.argsort(x[sel, 1])
        y = x[sel, 1][o] * 100
        ww = w[sel][o]
        uinf = 8.2 * np.sqrt(1.4)  # ADflow units: p_inf = rho_inf = 1
        p = 0.4 * (ww[:, 4] - 0.5 * ww[:, 0] * (ww[:, 1] ** 2 + ww[:, 2] ** 2))
        T = p / ww[:, 0]
        u = ww[:, 1] / uinf
        M = ww[:, 1] / np.sqrt(1.4 * T)
        main = m in ("SA", "SA-Edwards")
        for ax, v in zip(axs[:3], [u, T, M]):
            ax.plot(v, y, label=LABEL.get(m, m), lw=2.2 if main else 1.0, alpha=1.0 if main else 0.55,
                    zorder=3 if main else 2, **STYLE[m])
    for ax, key, xl in zip(axs[:3], ["U_over_Uinf", "T_over_Tinf", "M"], [r"$u/U_\infty$", r"$T/T_\infty$", "M"]):
        ax.plot(prof[key], prof["Y_cm"], label="Kussoy & Horstman survey (exp.)", zorder=4, **EXPSTYLE)
        ax.set_xlabel(xl)
        ax.set_ylabel("y [cm]")
        ax.set_ylim(0, 6)
        ax.grid(alpha=0.3)
    # Wall quantities and theta as bars
    ax = axs[3]
    c = e["conditions"]
    names = list(bl)
    quant = [("theta_cm", c["inflow_theta_cm"], r"$\theta$ [cm]"), ("Cf", c["inflow_Cf"], r"$C_f$"), ("qw_W_cm2", c["inflow_qw_W_cm2"], r"$q_w$ [W/cm$^2$]")]
    xpos = np.arange(len(quant))
    width = 0.8 / (len(names) + 1)
    for i, m in enumerate(names):
        r = bl[m]["x=-0.11"]
        vals = [r[k] / ref for k, ref, _ in quant]
        main = m in ("SA", "SA-Edwards")
        ax.bar(xpos + (i - len(names) / 2) * width, vals, width, color=STYLE[m]["color"], alpha=0.95 if main else 0.45,
               edgecolor="k" if main else "none", lw=0.8)
    ax.axhline(1.0, color="k", lw=1)
    ax.text(0.02, 0.97, "1.0 = experiment", transform=ax.transAxes, fontsize=8, va="top")
    ax.set_xticks(xpos)
    ax.set_xticklabels([q[2] for q in quant])
    ax.set_ylabel("ADflow / experiment")
    ax.set_title("ratio to the measured values at x = -11 cm\n" r"($\theta$ 0.094 cm*, $C_f$ 0.98e-3, $q_w$ 0.99 W/cm$^2$)", fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.03))
    fig.suptitle("Kussoy & Horstman Mach 8.2 flat plate, incoming boundary layer at x = -11 cm: ADflow (medium grid) vs experiment")
    fig.text(0.005, 0.005, "* theta matches by construction for SA: the plate leading edge is placed so its theta = 0.094 cm",
             ha="left", fontsize=7.5, color="0.35")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    savefig(fig, "plate_vs_experiment.png")


def plot_fin():
    res = load("fin")
    e = EXP["kussoy_fin15"]
    panels = [
        ("plate_p_X18.2", "Z_cm", "P_over_Pinf", "plate_p_X18.2_additional", "Plate pressure, X = 18.19 cm", "Z [cm]", r"$P/P_\infty$"),
        ("plate_q_X16.45", "Z_cm", "Q_over_Qinf", "plate_q_X16.45_additional", "Plate heat transfer, X = 16.45 cm", "Z [cm]", r"$Q_w/Q_{w,\infty}$"),
        ("plate_cf_X15.5", "Z_cm", "Cf", None, "Plate skin friction, X = 15.5 cm", "Z [cm]", r"$C_f$"),
        ("fin_p_X17.72", "Y_cm", "P_over_Pinf", None, "Fin pressure, S = 18.34 cm (X = 17.72 cm)", "Y [cm]", r"$P/P_\infty$"),
        ("fin_q_X16.13", "Y_cm", "Q_over_Qinf", "fin_q_SB_S19.62", "Fin heat transfer, S = 16.70 cm (X = 16.13 cm)", "Y [cm]", r"$Q_w/Q_{w,\infty}$"),
    ]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, (key, xk, yk, key2, title, xl, yl) in zip(axs.flat, panels):
        for m, d in res.items():
            ax.plot(d[key][xk], d[key]["value"], label=m, **STYLE[m])
        if key == "plate_cf_X15.5":
            ax.errorbar(e[key][xk], e[key][yk], yerr=e[key]["Cf_uncertainty"], capsize=3, label="Kussoy & Horstman (exp.)", **EXPSTYLE)
        else:
            ax.plot(e[key][xk], e[key][yk], label="Kussoy & Horstman (exp.)", **EXPSTYLE)
        if key2:
            lab = "exp., Schmidt-Boelter gauges, S = 19.62 cm" if key2.startswith("fin_q_SB") else "exp., additional data"
            ax.plot(e[key2][xk], e[key2][yk], label=lab, **EXPSTYLE2)
        ax.set_xlim(0, 15 if xk == "Z_cm" else 12)
        ax.set_ylim(bottom=0)
        ax.set_title(title)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axs.flat[-1].axis("off")
    c = e["conditions"]
    txt = [
        "Incoming boundary layer (x = -11 cm)",
        f"exp.: theta = {c['inflow_theta_cm']} cm, delta = {c['inflow_delta_cm']} cm",
    ]
    for m, d in res.items():
        b = d.get("incoming_BL", {})
        txt.append(f"{m}: theta = {b.get('theta_cm', np.nan):.3f} cm, d99 = {b.get('delta99_cm', np.nan):.2f} cm")
    txt += ["", "Experimental data: KUSSOY3.DAT (NASA TM-103838) as listed", "in Settles & Dodson, NASA CR-177638 (1994)."]
    axs.flat[-1].text(0.0, 0.9, "\n".join(txt), va="top", family="monospace", fontsize=9)
    fig.suptitle("Kussoy & Horstman Mach 8.2, 15 deg sharp fin: ADflow vs experiment")
    fig.tight_layout()
    savefig(fig, "fin_vs_experiment.png")


def plot_flare():
    res = load("flare")
    e = EXP["wideman_flare"]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, ph in zip(axs[0], [0, 90, 180]):
        for m, d in res.items():
            ax.plot(d[f"phi{ph}"]["x_cm"], d[f"phi{ph}"]["Cfx"], label=m, **STYLE[m])
        ax.plot(e[f"cfx_phi{ph}"]["x_cm"], e[f"cfx_phi{ph}"]["Cfx"], label="Wideman et al. (LISF)", **EXPSTYLE)
        if ph == 90:
            ax.plot(e["cfz_phi90_resolved"]["x_cm"], e["cfz_phi90_resolved"]["Cfx"], label="Wideman et al. (LISF, resolved)", **EXPSTYLE2)
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_xlim(-5, 20)
        ax.set_ylim(-0.001, 0.0045)
        ax.set_title(f"Streamwise skin friction, phi = {ph} deg")
        ax.set_xlabel("x [cm]")
        ax.set_ylabel(r"$C_{f,x}$")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    for ax, ph in zip(axs[1], [0, 90, 180]):
        for m, d in res.items():
            ax.plot(d[f"phi{ph}"]["x_cm"], d[f"phi{ph}"]["P_over_Pinf"], label=m, **STYLE[m])
        ax.plot(e["flare_pressure"][f"phi{ph}"]["x_cm"], e["flare_pressure"][f"phi{ph}"]["P_over_Pinf"], label="Wideman et al. (flare taps)", **EXPSTYLE)
        ax.plot(e["side_pressure"][f"phi{ph}"]["x_cm"], e["side_pressure"][f"phi{ph}"]["P_over_Pinf"], label="Wideman et al. (side taps)", **EXPSTYLE2)
        ax.set_xlim(-5, 25)
        ax.set_title(f"Surface pressure, phi = {ph} deg")
        ax.set_xlabel("x [cm]")
        ax.set_ylabel(r"$P/P_\infty$")
        ax.grid(alpha=0.3)
    axs[1, 0].legend(fontsize=7)
    fig.suptitle("Wideman et al. Mach 2.89 cylinder / offset 20 deg flare: ADflow vs experiment (NASA TM-108824 tables)")
    fig.tight_layout()
    savefig(fig, "flare_vs_experiment.png")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("plate", "all") and load("plate", "bl.json"):
        plot_plate()
    if which in ("fin", "all") and load("fin"):
        plot_fin()
    if which in ("flare", "all") and load("flare"):
        plot_flare()
