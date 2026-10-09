"""Boundary-layer profiles at x ~ -0.11 m: ADflow (cell centres) vs SU2 (nodes)."""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare import read_vtk, U, rho, Tinf, mu, run_dir, FIGURES
nuinf = mu / rho
level = sys.argv[1] if len(sys.argv) > 1 else "medium"  # needs results/adflow/plate/<level>/SA/cells.npz (dump_restart.py)
def profiles(model):
    d = np.load(os.path.join(run_dir("adflow", "plate", level, model), "cells.npz"))
    x, w = d["x"], d["w"]
    xu = np.unique(np.round(x[:, 0], 8)); xa = xu[np.argmin(np.abs(xu + 0.11))]
    m = (np.abs(x[:, 0] - xa) < 1e-7) & (x[:, 2] < x[:, 2].min() + 1e-9)
    o = np.argsort(x[m, 1]); yA = x[m, 1][o]; wA = w[m][o]
    winf1 = 8.2 * np.sqrt(1.4)  # U_inf in ADflow units (p_inf/rho_inf = 1)
    pA = 0.4 * (wA[:, 4] - 0.5 * wA[:, 0] * (wA[:, 1] ** 2 + wA[:, 2] ** 2))
    TA = pA / wA[:, 0]
    uA = wA[:, 1] / winf1
    nutA = wA[:, 5] / (wA[-1, 5] / 1.342)
    snap = sorted(glob.glob(os.path.join(run_dir("su2", "plate", level, model), "flow*_*.vtk")), key=os.path.getmtime)[-1]
    P, a = read_vtk(snap)
    xs = np.unique(np.round(P[:, 0], 9)); xb = xs[np.argmin(np.abs(xs - xa))]
    mm = np.abs(P[:, 0] - xb) < 1e-9; oo = np.argsort(P[mm, 1])
    yS = P[mm, 1][oo]; uS = (a["Momentum"][mm, 0] / a["Density"][mm])[oo] / U
    TS = a["Temperature"][mm][oo] / Tinf; nutS = a["Nu_Tilde"][mm][oo] / nuinf
    return xa, xb, (yA, uA, TA, nutA), (yS, uS, TS, nutS)
COL = {"SA": "k", "SA-Edwards": "tab:blue"}
fig, ax = plt.subplots(1, 3, figsize=(14, 4.9))
for model in ("SA", "SA-Edwards"):
    xa, xb, Ad, Sd = profiles(model)
    for k in range(3):
        ax[k].semilogx(Ad[0], Ad[k + 1], "-", color=COL[model], lw=1.6, label=f"{model}, ADflow")
        ax[k].semilogx(Sd[0], Sd[k + 1], "--", color=COL[model], lw=1.4, label=f"{model}, SU2")
    if model == "SA":
        yA, uA, TA, nutA = Ad; yS, uS, TS, nutS = Sd
for k, lab in enumerate((r"$u/U_\infty$", r"$T/T_\infty$", r"$\tilde\nu/\nu_\infty$")):
    ax[k].set_xlabel("wall distance y [m]"); ax[k].set_ylabel(lab); ax[k].grid(alpha=.3); ax[k].set_xlim(1e-6, 0.1)
ax[0].legend(fontsize=8.5, frameon=False, title="colour = model, line style = code", title_fontsize=8.5)
fig.suptitle(f"Mach 8.2 flat plate, boundary-layer profiles at x = {xa:.3f} m (ADflow cell centres; SU2 nodes at x = {xb:.3f} m), {level} grid")
fig.tight_layout(); ARCH = os.path.join(FIGURES, "..", "..", "archive", "figures_superseded"); os.makedirs(ARCH, exist_ok=True)   # superseded by the SU2 overlays in results/figures/by_model/plate_SA-Edwards.png
fig.savefig(os.path.join(ARCH, "plate_profiles_ADflow_vs_SU2" + ("" if level == "medium" else "_" + level) + ".png"), dpi=120)
for yy in [1e-4, 1e-3, 5e-3, 1e-2, 2e-2]:
    print(f"y={yy:.0e}: u/U A={np.interp(yy,yA,uA):.4f} S={np.interp(yy,yS,uS):.4f} | T/Tinf A={np.interp(yy,yA,TA):.3f} S={np.interp(yy,yS,TS):.3f} | nut/nu A={np.interp(yy,yA,nutA):.1f} S={np.interp(yy,yS,nutS):.1f}")
print("max nut/nu: ADflow %.1f SU2 %.1f" % (nutA.max(), nutS.max()))
