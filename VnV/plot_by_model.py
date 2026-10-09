"""One comparison figure per SA modification (each against baseline SA and the experiment), plus one bar
chart with all the models together.

Usage: MPLBACKEND=Agg python plot_by_model.py [plate|fin|flare|bars|all]
Writes results/figures/by_model/<case>_<model>.png and results/figures/plate_summary_bars.png.
Same data as plot_results.py (results/adflow/<case>/medium/<model>/, data/experiments_nasa.json).
"""
import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_results import EXP, STYLE, LABEL, EXPSTYLE, EXPSTYLE2, load, run_dir, FIGURES

OUT = os.path.join(FIGURES, "by_model")
DESC = {  # what each modification changes (README; all variants are SA-noft2)
    "SA-Edwards": "Edwards & Chandra (AIAA 94-2275): strain-based production, Shat = S (1/chi + fv1), "
                  "r = tanh(nuTilde/(Shat k^2 d^2))/tanh(1), no ft2",
    "SA-comp": "Spalart compressibility correction (AIAA 2000-2306): extra destruction "
               "-C5 nuTilde^2 (du_i/dx_j)^2 / a^2, C5 = 3.5",
    "SA-Catris": "Catris-Aupoix density correction (2000): diffuses sqrt(rho) nuTilde, laminar diffusion of mu, "
                 "original non-conservative left-hand side",
    "SA-Catris-cons": "Catris-Aupoix density correction, conservative left-hand side (TMR SA-noft2-CatrisCons)",
    "SA-Catris1999": "Catris (1999) form of the density-gradient diffusion terms",
    "SA-Seror": "Seror et al. (2005) form of the density-gradient diffusion terms",
    "SA-cons": "conservative diffusion form (TMR / SU2): adds the density-gradient cross terms",
    "SA-Edwards-Catris": "SA-Edwards combined with the Catris-Aupoix density correction",
}
BASE = dict(color="#2ca02c", ls="-", lw=2.0)                             # baseline SA: green, solid
MOD = dict(color="#d62728", ls="--", lw=2.4)                             # the modification: red, dashed
EXPA = dict(marker="o", mfc="none", mec="k", color="k", mew=1.4, ls="none", ms=5.5)   # main data set
EXPB = dict(marker="s", mfc="none", mec="0.45", color="0.45", mew=1.2, ls="none", ms=4.5)   # second data set: grey
UINF = 8.2 * np.sqrt(1.4)          # ADflow units (p_inf = rho_inf = 1)


def plate_profile(m):
    d = np.load(os.path.join(run_dir("adflow", "plate", "medium", m), "cells.npz"))
    x, w = d["x"], d["w"]
    xu = np.unique(np.round(x[:, 0], 8))
    xa = xu[np.argmin(np.abs(xu + 0.11))]
    sel = (np.abs(x[:, 0] - xa) < 1e-7) & (x[:, 2] < x[:, 2].min() + 1e-9)
    o = np.argsort(x[sel, 1])
    y = x[sel, 1][o] * 100                                   # cm, wall at y = 0
    ww = w[sel][o]
    p = 0.4 * (ww[:, 4] - 0.5 * ww[:, 0] * (ww[:, 1] ** 2 + ww[:, 2] ** 2))
    T = p / ww[:, 0]
    return dict(y=y, u=ww[:, 1] / UINF, T=T, M=ww[:, 1] / np.sqrt(1.4 * T), nut=ww[:, 5] / (ww[-1, 5] / 1.342))


def plot_plate():
    e = EXP["kussoy_fin15"]; prof = e["upstream_BL_profile"]; c = e["conditions"]
    bl = load("plate", "bl.json")
    P = {m: plate_profile(m) for m in bl}
    exp_int = {"theta_cm": c["inflow_theta_cm"], "Cf": c["inflow_Cf"], "qw_W_cm2": c["inflow_qw_W_cm2"]}
    for m in bl:
        if m == "SA":
            continue
        st = MOD
        fig, axs = plt.subplots(2, 4, figsize=(17, 8.4))
        for row, logy in enumerate((False, True)):
            for k, (key, ek, xl) in enumerate((("u", "U_over_Uinf", r"$u/U_\infty$"), ("T", "T_over_Tinf", r"$T/T_\infty$"),
                                               ("M", "M", "Mach"), ("nut", None, r"$\tilde\nu/\nu_\infty$"))):
                if row == 0 and key == "nut":
                    continue
                ax = axs[row, k]
                for mm, s, lab in (("SA", BASE, "SA (baseline)"), (m, st, LABEL.get(m, m))):
                    yy = np.maximum(P[mm]["y"], 1e-6) if logy else P[mm]["y"]
                    ax.plot(P[mm][key], yy, label=lab, **s)
                if ek:
                    ye = np.array(prof["Y_cm"]); ve = np.array(prof[ek]); keep = ye > 0 if logy else ye >= 0
                    ax.plot(ve[keep], ye[keep], label="Kussoy & Horstman survey", zorder=4, **EXPA)
                ax.set_xlabel(xl); ax.grid(alpha=0.3)
                if logy:
                    ax.set_yscale("log"); ax.set_ylim(1e-4, 6); ax.set_ylabel("wall distance y [cm] (log)")
                else:
                    ax.set_ylim(0, 6); ax.set_ylabel("y [cm]")
        for axx in list(axs[0, :3]) + list(axs[1, :]):
            axx.legend(fontsize=7.5, frameon=True, framealpha=0.85, loc="best")
        axs[0, 0].set_title("linear height (whole layer)", loc="left", fontsize=9, color="0.35")
        axs[1, 0].set_title("log wall distance (near wall; T peaks then falls to the 300 K wall)", loc="left", fontsize=9, color="0.35")
        # integral / wall quantities table
        ax = axs[0, 3]; ax.axis("off")
        rows = [("quantity", "exp.", "SA", "model", "model\nvs SA", "model\nvs exp.")]
        for k, lab, f in (("theta_cm", "theta, cm", "{:.4f}"), ("delta99_cm", "d99, cm", "{:.3f}"),
                          ("Cf", "Cf x 1e3", "{:.3f}"), ("qw_W_cm2", "q_w, W/cm2", "{:.3f}")):
            a, b = bl["SA"]["x=-0.11"][k], bl[m]["x=-0.11"][k]
            ev = exp_int.get(k)
            if k == "Cf":
                a, b, ev = a * 1e3, b * 1e3, ev * 1e3
            rows.append((lab, f.format(ev) if ev else "-", f.format(a), f.format(b), f"{100 * (b / a - 1):+.1f}%",
                         f"{100 * (b / ev - 1):+.1f}%" if ev else "-"))
        ax.set_position([0.765, 0.55, 0.228, 0.33])
        tb = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center",
                      colWidths=[0.24, 0.15, 0.15, 0.15, 0.155, 0.155])
        tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1.0, 2.0)
        for (r, cc), cell in tb.get_celld().items():
            if r == 0:
                cell.set_text_props(weight="bold"); cell.set_height(cell.get_height() * 1.5)
        ax.set_title(f"x = -11 cm; model = {m}\n(theta matches for SA by construction)", fontsize=9)
        fig.suptitle(f"Mach 8.2 flat plate, incoming boundary layer at x = -11 cm:  {LABEL.get(m, m)}  vs  SA  vs  experiment\n"
                     f"{DESC.get(m, '')}", fontsize=11)
        fig.subplots_adjust(left=0.05, right=0.99, bottom=0.07, top=0.88, wspace=0.28, hspace=0.3)
        ax.set_position([0.765, 0.55, 0.228, 0.33])
        os.makedirs(OUT, exist_ok=True); fig.savefig(os.path.join(OUT, f"plate_{m}.png"), dpi=120); plt.close(fig)
        print("wrote", os.path.join(OUT, f"plate_{m}.png"))


def plot_bars():
    e = EXP["kussoy_fin15"]; c = e["conditions"]; bl = load("plate", "bl.json"); names = list(bl)
    quant = [("theta_cm", c["inflow_theta_cm"], r"$\theta$ (matched for SA)"), ("Cf", c["inflow_Cf"], r"$C_f$"),
             ("qw_W_cm2", c["inflow_qw_W_cm2"], r"$q_w$")]
    fig, ax = plt.subplots(figsize=(13, 5.2))
    xpos = np.arange(len(quant)); width = 0.85 / len(names)
    for i, m in enumerate(names):
        vals = [bl[m]["x=-0.11"][k] / ref for k, ref, _ in quant]
        xs = xpos + (i - (len(names) - 1) / 2) * width
        ax.bar(xs, vals, width, color=BASE["color"] if m == "SA" else STYLE[m]["color"], edgecolor="k", lw=0.5,
               label=LABEL.get(m, m))
        for xx, v in zip(xs, vals):
            ax.text(xx, v + 0.01, f"{v:.2f}", ha="center", va="bottom", fontsize=6.5, rotation=90)
    ax.axhline(1.0, color="k", lw=1)
    ax.set_xticks(xpos); ax.set_xticklabels([q[2] for q in quant], fontsize=11)
    ax.set_ylabel("ADflow / experiment  (1 = experiment)"); ax.set_ylim(0, 1.35); ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8.5, ncol=3, frameon=False, loc="upper left")
    ax.set_title(f"Mach 8.2 flat plate at x = -11 cm, all models (medium grid). Experiment: theta {c['inflow_theta_cm']} cm, "
                 f"Cf {c['inflow_Cf']:.2e}, q_w {c['inflow_qw_W_cm2']} W/cm$^2$", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIGURES, "plate_summary_bars.png"), dpi=130); plt.close(fig)
    print("wrote", os.path.join(FIGURES, "plate_summary_bars.png"))


FIN_PANELS = [
    ("plate_p_X18.2", "Z_cm", "P_over_Pinf", "plate_p_X18.2_additional", "Plate pressure, X = 18.19 cm", "Z [cm]", r"$P/P_\infty$"),
    ("plate_q_X16.45", "Z_cm", "Q_over_Qinf", "plate_q_X16.45_additional", "Plate heat transfer, X = 16.45 cm", "Z [cm]", r"$Q_w/Q_{w,\infty}$"),
    ("plate_cf_X15.5", "Z_cm", "Cf", None, "Plate skin friction, X = 15.5 cm", "Z [cm]", r"$C_f$"),
    ("fin_p_X17.72", "Y_cm", "P_over_Pinf", None, "Fin pressure, S = 18.34 cm", "Y [cm]", r"$P/P_\infty$"),
    ("fin_q_X16.13", "Y_cm", "Q_over_Qinf", "fin_q_SB_S19.62", "Fin heat transfer, S = 16.70 cm", "Y [cm]", r"$Q_w/Q_{w,\infty}$"),
]


def plot_fin():
    res = load("fin"); e = EXP["kussoy_fin15"]
    for m in res:
        if m == "SA":
            continue
        st = MOD
        fig, axs = plt.subplots(2, 3, figsize=(15, 8.6))
        for ax, (key, xk, yk, key2, title, xl, yl) in zip(axs.flat, FIN_PANELS):
            ax.plot(res["SA"][key][xk], res["SA"][key]["value"], label="SA (baseline)", **BASE)
            ax.plot(res[m][key][xk], res[m][key]["value"], label=LABEL.get(m, m), **st)
            if key == "plate_cf_X15.5":
                ax.errorbar(e[key][xk], e[key][yk], yerr=e[key]["Cf_uncertainty"], capsize=3, label="Kussoy & Horstman", **EXPA)
            else:
                ax.plot(e[key][xk], e[key][yk], label="Kussoy & Horstman", **EXPA)
            if key2:
                ax.plot(e[key2][xk], e[key2][yk], label="exp., S-B gauges S = 19.62 cm" if key2.startswith("fin_q_SB") else "exp., second set", **EXPB)
            ax.set_xlim(0, 15 if xk == "Z_cm" else 12); ax.set_ylim(bottom=0)
            ax.set_title(title, fontsize=10); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=0.3)
            ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, loc="best")
        axs.flat[-1].axis("off")
        fig.suptitle(f"Kussoy & Horstman Mach 8.2, 15 deg sharp fin:  {LABEL.get(m, m)}  vs  SA  vs  experiment\n{DESC.get(m, '')}", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        os.makedirs(OUT, exist_ok=True); fig.savefig(os.path.join(OUT, f"fin_{m}.png"), dpi=120); plt.close(fig)
        print("wrote", os.path.join(OUT, f"fin_{m}.png"))


def plot_flare():
    res = load("flare"); e = EXP["wideman_flare"]
    for m in res:
        if m == "SA":
            continue
        st = MOD
        fig, axs = plt.subplots(2, 3, figsize=(15, 8.6))
        for ax, ph in zip(axs[0], [0, 90, 180]):
            ax.plot(res["SA"][f"phi{ph}"]["x_cm"], res["SA"][f"phi{ph}"]["Cfx"], label="SA (baseline)", **BASE)
            ax.plot(res[m][f"phi{ph}"]["x_cm"], res[m][f"phi{ph}"]["Cfx"], label=LABEL.get(m, m), **st)
            ax.plot(e[f"cfx_phi{ph}"]["x_cm"], e[f"cfx_phi{ph}"]["Cfx"], label="Wideman et al. (LISF)", **EXPA)
            if ph == 90:
                ax.plot(e["cfz_phi90_resolved"]["x_cm"], e["cfz_phi90_resolved"]["Cfx"], label="Wideman et al. (LISF, resolved)", **EXPB)
            ax.axhline(0, color="gray", lw=0.5); ax.set_xlim(-5, 20); ax.set_ylim(-0.001, 0.0045)
            ax.set_title(f"Streamwise skin friction, phi = {ph} deg", fontsize=10); ax.set_xlabel("x [cm]"); ax.set_ylabel(r"$C_{f,x}$")
            ax.grid(alpha=0.3); ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, loc="best")
        for ax, ph in zip(axs[1], [0, 90, 180]):
            ax.plot(res["SA"][f"phi{ph}"]["x_cm"], res["SA"][f"phi{ph}"]["P_over_Pinf"], label="SA (baseline)", **BASE)
            ax.plot(res[m][f"phi{ph}"]["x_cm"], res[m][f"phi{ph}"]["P_over_Pinf"], label=LABEL.get(m, m), **st)
            ax.plot(e["flare_pressure"][f"phi{ph}"]["x_cm"], e["flare_pressure"][f"phi{ph}"]["P_over_Pinf"], label="Wideman et al. (flare taps)", **EXPA)
            ax.plot(e["side_pressure"][f"phi{ph}"]["x_cm"], e["side_pressure"][f"phi{ph}"]["P_over_Pinf"], label="Wideman et al. (side taps)", **EXPB)
            ax.set_xlim(-5, 25); ax.set_title(f"Surface pressure, phi = {ph} deg", fontsize=10); ax.set_xlabel("x [cm]")
            ax.set_ylabel(r"$P/P_\infty$"); ax.grid(alpha=0.3); ax.legend(fontsize=7.5, frameon=True, framealpha=0.9, loc="best")
        fig.suptitle(f"Wideman et al. Mach 2.89 cylinder / offset 20 deg flare:  {LABEL.get(m, m)}  vs  SA  vs  experiment\n{DESC.get(m, '')}", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        os.makedirs(OUT, exist_ok=True); fig.savefig(os.path.join(OUT, f"flare_{m}.png"), dpi=120); plt.close(fig)
        print("wrote", os.path.join(OUT, f"flare_{m}.png"))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("plate", "all"):
        plot_plate()
    if which in ("bars", "all"):
        plot_bars()
    if which in ("fin", "all"):
        plot_fin()
    if which in ("flare", "all"):
        plot_flare()
