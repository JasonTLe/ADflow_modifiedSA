"""Compare ADflow and SU2 on the 2-D Mach 8.2 cold-wall flat plate.

Usage: MPLBACKEND=Agg python compare.py
For each grid level (coarse, medium, fine) and model it reads
  ADflow: results/adflow/plate/<level>/<model>/*_surf.plt, bl.json
  SU2:    results/su2/plate/<level>/<model>/flow*_*.vtk
and writes results/su2_comparison.json, results/su2_comparison.md and
results/figures/plate_ADflow_vs_SU2.png. Wall quantities are compared
as functions of x; heat flux uses the same reference for both codes
(q = St * rho_inf U_inf cp (T0 - Tw), cp = gamma R / (gamma - 1))."""

import os
import sys
import glob
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from common import read_tecplot_surface, face_centers, run_dir, RESULTS, FIGURES  # noqa: E402

M, Tinf, Re, R, gam, Tw = 8.2, 76.89, 5.32e6, 287.05, 1.4, 300.0
mu = 1.716e-5 * (Tinf / 273.15) ** 1.5 * (273.15 + 110.55) / (Tinf + 110.55)
a = np.sqrt(gam * R * Tinf)
U = M * a
rho = Re * mu / U
cp = gam * R / (gam - 1)
T0 = Tinf * (1 + 0.5 * (gam - 1) * M**2)
qref = rho * U * cp * (T0 - Tw)  # W/m^2 per unit Stanton number
qdyn = 0.5 * rho * U**2

LEVELS = ["coarse", "medium", "fine"]
MODELS = ["SA", "SA-Edwards"]
STATIONS = [-0.5, -0.11, 0.0]


def adflow(level, model):
    d = run_dir("adflow", "plate", level, model)
    f = sorted(glob.glob(os.path.join(d, "*_surf.plt")), key=os.path.getmtime)
    if not f:
        return None
    z = read_tecplot_surface(f[-1])["plate"]
    X, cf = face_centers(z, "SkinFrictionMagnitude")
    _, st = face_centers(z, "StantonNumber")
    m = X[:, 2] < X[:, 2].min() + 1e-9  # one element row in z
    o = np.argsort(X[m, 0])
    bl = json.load(open(os.path.join(d, "bl.json")))
    return {"x": X[m, 0][o], "Cf": cf[m][o], "q": st[m][o] * qref, "bl": bl}


def read_vtk(fileName):
    """Minimal reader for SU2's legacy ASCII VTK volume output."""
    L = open(fileName).read().split("\n")
    i = next(k for k, l in enumerate(L) if l.startswith("POINTS"))
    n = int(L[i].split()[1])

    def block(start, count):
        v, j = [], start
        while len(v) < count:
            v += L[j].split()
            j += 1
        return np.array(v[:count], float)

    P = block(i + 1, 3 * n).reshape(n, 3)
    arr = {}
    for k, l in enumerate(L):
        if l.startswith("SCALARS"):
            arr[l.split()[1]] = block(k + 2, n)
        elif l.startswith("VECTORS"):
            arr[l.split()[1]] = block(k + 1, 3 * n).reshape(n, 3)
    return P, arr


def su2_fields(fileName):
    P, a = read_vtk(fileName)
    x, y = P[:, 0], P[:, 1]
    wall = (np.abs(y) < 1e-12) & (x > -1.64 + 1e-9)
    o = np.argsort(x[wall])
    res = {
        "x": x[wall][o],
        "Cf": np.abs(a["Skin_Friction_Coefficient"][wall, 0][o]),
        "q": np.abs(a["Heat_Flux"][wall][o]),
    }
    rho_, u = a["Density"], a["Momentum"][:, 0] / a["Density"]
    ux = np.unique(np.round(x, 9))
    th = []
    for xl in ux:
        mm = np.abs(x - xl) < 1e-9
        oo = np.argsort(y[mm])
        yy, rr, uu = y[mm][oo], rho_[mm][oo], u[mm][oo]
        i99 = np.argmax(uu >= 0.99 * U)
        d99 = np.interp(0.99 * U, uu[i99 - 1 : i99 + 1], yy[i99 - 1 : i99 + 1]) if i99 > 0 else np.nan
        k = yy <= 1.5 * d99
        f2 = rr * uu / (rho * U) * (1 - uu / U)
        th.append(np.trapz(f2[k], yy[k]) if k.sum() > 1 else np.nan)
    res["theta_x"] = ux
    res["theta"] = np.array(th)
    R = a["Residual_Density"]
    away = x > -1.60
    res["log_res_rho_all"] = float(np.log10(np.sqrt(np.mean(R**2))))
    res["log_res_rho_awayLE"] = float(np.log10(np.sqrt(np.mean(R[away] ** 2))))
    return res


def su2(level, model):
    d = run_dir("su2", "plate", level, model)
    # latest snapshots (stage 2 writes flow_*.vtk, a continuation run flow3_*.vtk)
    snaps = sorted(glob.glob(os.path.join(d, "flow*_*.vtk")), key=os.path.getmtime)
    if not os.path.exists(os.path.join(d, "stage2.log")) or not snaps:
        return None
    res = su2_fields(snaps[-1])
    res["snapshot"] = os.path.basename(snaps[-1])
    if len(snaps) > 1:
        prev = su2_fields(snaps[-2])
        for q in ["Cf", "q"]:
            res[f"d{q}_last_snapshots_pct"] = float(
                100 * np.max(np.abs(res[q] - prev[q])[res["x"] > -1.5]) / np.max(res[q][res["x"] > -1.5])
            )
    return res


out = {}
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5))
colors = {"coarse": "tab:green", "medium": "k", "fine": "tab:blue", "finer": "tab:purple"}
for name in LEVELS:
    for model in MODELS:
        A, S = adflow(name, model), su2(name, model)
        key = f"{model}/{name}"
        out[key] = {}
        ls = "-" if model == "SA" else "--"
        col = colors[name]
        row = 0 if model == "SA" else 1
        if A is not None:
            axs[row, 0].plot(A["x"], A["Cf"] * 1e3, ls, color=col, lw=1.2, label=f"ADflow {name}")
            axs[row, 1].plot(A["x"], A["q"] / 1e4, ls, color=col, lw=1.2, label=f"ADflow {name}")
        if S is not None:
            axs[row, 0].plot(S["x"], S["Cf"] * 1e3, ":", color=col, lw=2, label=f"SU2 {name}")
            axs[row, 1].plot(S["x"], S["q"] / 1e4, ":", color=col, lw=2, label=f"SU2 {name}")
        for xs in STATIONS:
            e = {}
            if A is not None:
                e["ADflow_Cf"] = float(np.interp(xs, A["x"], A["Cf"]))
                e["ADflow_qw_W_cm2"] = float(np.interp(xs, A["x"], A["q"]) / 1e4)
                blk = f"x={xs:+.2f}"
                if blk in A["bl"]:
                    e["ADflow_theta_cm"] = A["bl"][blk]["theta_cm"]
                    e["ADflow_theta_x"] = A["bl"][blk]["x"]
            if S is not None:
                e["SU2_Cf"] = float(np.interp(xs, S["x"], S["Cf"]))
                e["SU2_qw_W_cm2"] = float(np.interp(xs, S["x"], S["q"]) / 1e4)
                xt = e.get("ADflow_theta_x", xs)
                e["SU2_theta_cm"] = float(np.interp(xt, S["theta_x"], S["theta"]) * 100)
                e["SU2_log_res_rho_all"] = S["log_res_rho_all"]
                e["SU2_log_res_rho_awayLE"] = S["log_res_rho_awayLE"]
                e["SU2_dCf_last_snapshots_pct"] = S.get("dCf_last_snapshots_pct")
                e["SU2_dq_last_snapshots_pct"] = S.get("dq_last_snapshots_pct")
            if "ADflow_Cf" in e and "SU2_Cf" in e:
                e["dCf_pct"] = 100 * (e["SU2_Cf"] / e["ADflow_Cf"] - 1)
                e["dq_pct"] = 100 * (e["SU2_qw_W_cm2"] / e["ADflow_qw_W_cm2"] - 1)
                if "ADflow_theta_cm" in e:
                    e["dtheta_pct"] = 100 * (e["SU2_theta_cm"] / e["ADflow_theta_cm"] - 1)
            out[key][f"x={xs:+.2f}"] = e
for row, model in enumerate(MODELS):
    axs[row, 0].set_ylabel(f"{model}: $C_f$ x 1e3")
    axs[row, 1].set_ylabel(f"{model}: $q_w$ [W/cm$^2$]")
    for c in range(2):
        axs[row, c].set_xlim(-1.64, 0.25)
        axs[row, c].set_xlabel("x [m] (plate leading edge at -1.64 m)")
        axs[row, c].grid(alpha=0.3)
    axs[row, 0].set_ylim(0.5, 2.0)
    axs[row, 1].set_ylim(0.5, 2.0)
    axs[row, 0].plot([-0.11], [1.0], "ko", mfc="none", label="Kussoy exp. (x = -11 cm)")
    axs[row, 1].plot([-0.11], [1.04], "ko", mfc="none", label="Kussoy exp. (x = -11 cm)")
axs[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("Mach 8.2 cold-wall flat plate: ADflow (solid/dashed) vs SU2 (dotted), same grids")
fig.tight_layout()
os.makedirs(FIGURES, exist_ok=True)
fig.savefig(os.path.join(FIGURES, "plate_ADflow_vs_SU2.png"), dpi=130)
json.dump(out, open(os.path.join(RESULTS, "su2_comparison.json"), "w"), indent=1)

# Markdown tables at x = -11 cm
X = "x=-0.11"
L = ["### ADflow vs SU2, Mach 8.2 flat plate, x = -11 cm", "",
     "| model | grid | Cf x1e3 ADflow / SU2 | q_w [W/cm^2] ADflow / SU2 | theta [cm] ADflow / SU2 | SU2 - ADflow (Cf / q_w / theta) |",
     "|---|---|---|---|---|---|"]
for model in MODELS:
    for name in LEVELS:
        e = out[f"{model}/{name}"].get(X, {})
        if "dCf_pct" not in e:
            continue
        L.append(f"| {model} | {name} | {e['ADflow_Cf']*1e3:.3f} / {e['SU2_Cf']*1e3:.3f} "
                 f"| {e['ADflow_qw_W_cm2']:.3f} / {e['SU2_qw_W_cm2']:.3f} "
                 f"| {e['ADflow_theta_cm']:.4f} / {e['SU2_theta_cm']:.4f} "
                 f"| {e['dCf_pct']:+.1f} / {e['dq_pct']:+.1f} / {e['dtheta_pct']:+.1f} % |")
L += ["", "### SA-Edwards relative to SA (x = -11 cm)", "",
      "| grid | Cf ADflow / SU2 | q_w ADflow / SU2 | theta ADflow / SU2 |", "|---|---|---|---|"]
for name in LEVELS:
    a, b = out[f"SA/{name}"].get(X, {}), out[f"SA-Edwards/{name}"].get(X, {})
    if "dCf_pct" not in a or "dCf_pct" not in b:
        continue
    r = lambda k: 100 * (b[k] / a[k] - 1)
    L.append(f"| {name} | {r('ADflow_Cf'):+.2f} / {r('SU2_Cf'):+.2f} % | {r('ADflow_qw_W_cm2'):+.2f} / {r('SU2_qw_W_cm2'):+.2f} % "
             f"| {r('ADflow_theta_cm'):+.1f} / {r('SU2_theta_cm'):+.1f} % |")
open(os.path.join(RESULTS, "su2_comparison.md"), "w").write("\n".join(L) + "\n")
for k, v in out.items():
    e = v.get("x=-0.11", {})
    if e:
        print(k, {kk: round(vv, 5) if isinstance(vv, float) else vv for kk, vv in e.items()})
