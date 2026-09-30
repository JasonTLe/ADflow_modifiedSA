"""Plot ADflow results of the validation cases against the experimental data.

Usage: MPLBACKEND=Agg python plot_results.py [fin|flare|all]
Reads results/adflow/<case>/medium/<model>/lines.json and data/experiments.json
and writes results/figures/<case>_vs_experiment.png."""

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
EXP = json.load(open(os.path.join(HERE, "data", "experiments.json")))
STYLE = {
    "SA": dict(color="k", ls="-"),
    "SA-Edwards": dict(color="tab:blue", ls="-"),
    "SA-comp": dict(color="tab:purple", ls="--"),
}


def load(case):
    out = {}
    for f in sorted(glob.glob(os.path.join(run_dir("adflow", case, "medium", "*"), "lines.json"))):
        d = json.load(open(f))
        name = os.path.basename(os.path.dirname(f))
        if name in STYLE:
            out[name] = d
    return out


def plot_fin():
    res = load("fin")
    e = EXP["kussoy_fin15"]
    panels = [
        ("plate_p_X18.2", "Z_cm", "P_over_Pinf", "Plate pressure, X = 18.2 cm", "Z [cm]", r"$P/P_\infty$"),
        ("plate_q_X16.45", "Z_cm", "Q_over_Qinf", "Plate heat transfer, X = 16.45 cm", "Z [cm]", r"$Q_w/Q_{w,\infty}$"),
        ("plate_cf_X15.5", "Z_cm", "Cf", "Plate skin friction, X = 15.5 cm", "Z [cm]", r"$C_f$"),
        ("fin_p_X17.72", "Y_cm", "P_over_Pinf", "Fin pressure, X = 17.72 cm", "Y [cm]", r"$P/P_\infty$"),
        ("fin_q_X16.13", "Y_cm", "Q_over_Qinf", "Fin heat transfer, X = 16.13 cm", "Y [cm]", r"$Q_w/Q_{w,\infty}$"),
    ]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, (key, xk, yk, title, xl, yl) in zip(axs.flat, panels):
        for m, d in res.items():
            ax.plot(d[key][xk], d[key]["value"], label=m, **STYLE[m])
        ax.plot(e[key][xk], e[key][yk], "o", mfc="none", color="k", label="Kussoy & Horstman (exp.)")
        ax.set_xlim(0, 15 if xk == "Z_cm" else 12)
        ax.set_ylim(bottom=0)
        ax.set_title(title)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
    axs.flat[0].legend(fontsize=8)
    axs.flat[-1].axis("off")
    txt = ["Incoming boundary layer (x = -11 cm)", "exp.: theta = 0.094 cm, delta = 3.7 cm"]
    for m, d in res.items():
        b = d.get("incoming_BL", {})
        txt.append(f"{m}: theta = {b.get('theta_cm', np.nan):.3f} cm, d99 = {b.get('delta99_cm', np.nan):.2f} cm")
    axs.flat[-1].text(0.0, 0.9, "\n".join(txt), va="top", family="monospace", fontsize=9)
    fig.suptitle("Kussoy & Horstman Mach 8.2, 15 deg sharp fin: ADflow vs experiment")
    fig.tight_layout()
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, "fin_vs_experiment.png"), dpi=130)


def plot_flare():
    res = load("flare")
    e = EXP["wideman_flare"]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    for ax, ph in zip(axs[0], [0, 90, 180]):
        for m, d in res.items():
            ax.plot(d[f"phi{ph}"]["x_cm"], d[f"phi{ph}"]["Cfx"], label=m, **STYLE[m])
        ax.plot(e[f"cfx_phi{ph}"]["x_cm"], e[f"cfx_phi{ph}"]["Cfx"], "o", mfc="none", color="k", label="Wideman et al. (exp.)")
        ax.axhline(0, color="gray", lw=0.5)
        ax.set_xlim(-5, 20)
        ax.set_ylim(-0.001, 0.0045)
        ax.set_title(f"Streamwise skin friction, phi = {ph} deg")
        ax.set_xlabel("x [cm]")
        ax.set_ylabel(r"$C_{f,x}$")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=8)
    for ax, ph in zip(axs[1], [0, 90, 180]):
        for m, d in res.items():
            ax.plot(d[f"phi{ph}"]["x_cm"], d[f"phi{ph}"]["P_over_Pinf"], label=m, **STYLE[m])
        ax.set_xlim(-5, 25)
        ax.set_title(f"Surface pressure, phi = {ph} deg")
        ax.set_xlabel("x [cm]")
        ax.set_ylabel(r"$P/P_\infty$")
        ax.grid(alpha=0.3)
    t = e["text_metrics"]
    axs[1, 0].axvline(t["phi0"]["p_plateau_x_cm"], color="gray", ls=":")
    axs[1, 0].axhline(t["phi0"]["p_rise_ratio"], color="gray", ls="--", lw=0.8)
    axs[1, 0].text(-4.5, t["phi0"]["p_rise_ratio"] + 0.03, "exp. rise ratio 2.59, plateau at x = 2.35 cm", fontsize=8)
    axs[1, 2].text(-4.5, 2.6, "exp. plateau at x = 9.25 cm (dotted)", fontsize=8)
    axs[1, 2].axvline(t["phi180"]["p_plateau_x_cm"], color="gray", ls=":")
    fig.suptitle("Wideman et al. Mach 2.89 cylinder / offset 20 deg flare: ADflow vs experiment")
    fig.tight_layout()
    os.makedirs(FIGURES, exist_ok=True)
    fig.savefig(os.path.join(FIGURES, "flare_vs_experiment.png"), dpi=130)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("fin", "all") and load("fin"):
        plot_fin()
    if which in ("flare", "all") and load("flare"):
        plot_flare()
