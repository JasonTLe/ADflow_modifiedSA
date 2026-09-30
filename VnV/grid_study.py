"""Grid-convergence study (Roache GCI) of the validation cases.

Usage: python grid_study.py  -> prints Markdown tables, writes results/grid_study.md
Levels: 3-D cases refined by r = sqrt(2) per direction (refc, medium, reff);
2-D plate by r = 2 (refc, medium, reff).

For three solutions f1 (fine), f2, f3 (coarse) with constant ratio r:
  p = ln(|(f3 - f2) / (f2 - f1)|) / ln(r)           observed order
  f_ext = f1 + (f1 - f2) / (r^p - 1)                  Richardson extrapolation
  GCI_fine = 1.25 |(f2 - f1) / f1| / (r^p - 1)        (Roache, Fs = 1.25)
If (f3 - f2) and (f2 - f1) have opposite signs the convergence is oscillatory;
then p is not defined and the bound max|f_i - f_j| / |f1| is reported.
p is limited to [0.5, 3] when computing GCI."""

import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from common import RESULTS, run_dir  # noqa: E402
R3D, R2D = np.sqrt(2.0), 2.0


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def gci(f, r):
    """f = [coarse, medium, fine]."""
    f3, f2, f1 = f
    e21, e32 = f1 - f2, f2 - f3
    out = {"coarse": f3, "medium": f2, "fine": f1}
    if e21 == 0:
        out.update(p=np.nan, ext=f1, gci_pct=0.0, note="identical")
        return out
    if e32 * e21 < 0:
        out.update(p=np.nan, ext=np.nan, gci_pct=100 * max(abs(f1 - f2), abs(f2 - f3), abs(f1 - f3)) / abs(f1), note="oscillatory")
        return out
    p = np.log(abs(e32 / e21)) / np.log(r)
    pl = float(np.clip(p, 0.5, 3.0))
    out.update(p=float(p), ext=f1 + e21 / (r**pl - 1), gci_pct=100 * 1.25 * abs(e21 / f1) / (r**pl - 1), note="monotonic")
    return out


def band(d, key, lo, hi, xk):
    x = np.array(d[key][xk])
    v = np.array(d[key]["value"])
    return float(np.mean(v[(x >= lo) & (x <= hi)]))


def fin_metrics(d):
    return {
        "plate P peak (X=18.2)": max(d["plate_p_X18.2"]["value"]),
        "plate P plateau (Z=4-6)": band(d, "plate_p_X18.2", 4, 6, "Z_cm"),
        "plate Q peak (X=16.45)": max(d["plate_q_X16.45"]["value"]),
        "plate Cf peak x1e3 (X=15.5)": 1e3 * max(d["plate_cf_X15.5"]["value"]),
        "fin P (Y=6-10)": band(d, "fin_p_X17.72", 6, 10, "Y_cm"),
        "fin Q (Y=6-10)": band(d, "fin_q_X16.13", 6, 10, "Y_cm"),
        "incoming theta [cm]": d["incoming_BL"]["theta_cm"],
    }


def flare_metrics(d):
    p0, p90, p180 = d["phi0"], d["phi90"], d["phi180"]
    x0 = np.array(p0["x_cm"])
    zc = p180["Cfx_zero_crossings_cm"]
    x180, c180 = np.array(p180["x_cm"]), np.array(p180["Cfx"])
    x90, c90 = np.array(p90["x_cm"]), np.array(p90["Cfx"])
    return {
        "Cf upstream x1e3": 1e3 * p0["Cfx_upstream_x-10cm"],
        "phi=0 P rise": float(np.max(np.array(p0["P_over_Pinf"])[(x0 > 0) & (x0 < 10)]) / p0["P_upstream_x-10cm"]),
        "phi=180 x_sep [cm]": zc[0] if zc else np.nan,
        "phi=180 x_att [cm]": zc[1] if len(zc) > 1 else np.nan,
        "phi=180 Cf min x1e3": 1e3 * c180[(x180 > 3) & (x180 < 9)].min(),
        "phi=90 Cf max x1e3": 1e3 * c90[(x90 > 5) & (x90 < 14)].max(),
        "delta99 x=-2cm [cm]": d["incoming_BL"]["delta99_cm"],
    }


def plate_metrics(d, runDir):
    """theta from bl.json (nearest cell line, x within 3 mm of -0.11 m; theta
    varies by < 0.01 % over that distance); Cf and q_w interpolated to
    exactly x = -0.11 m from the surface solution."""
    import glob
    from common import read_tecplot_surface, face_centers

    r = d["x=-0.11"]
    z = read_tecplot_surface(sorted(glob.glob(os.path.join(runDir, "*_surf.plt")), key=os.path.getmtime)[-1])["plate"]
    X, cf = face_centers(z, "SkinFrictionMagnitude")
    _, st = face_centers(z, "StantonNumber")
    m = X[:, 2] < X[:, 2].min() + 1e-9
    o = np.argsort(X[m, 0])
    x = X[m, 0][o]
    q_per_st = r["qw_W_cm2"] / np.interp(r["x"], x, st[m][o])  # same conversion as run_flatplate.py
    return {
        "theta [cm]": r["theta_cm"],
        "Cf x1e3": 1e3 * float(np.interp(-0.11, x, cf[m][o])),
        "q_w [W/cm2]": float(np.interp(-0.11, x, st[m][o]) * q_per_st),
    }


def table(title, rows, r, labels=("coarse", "medium", "fine")):
    lines = [f"### {title}\n", f"| Quantity | {labels[0]} | {labels[1]} | {labels[2]} | observed p | extrapolated | GCI_fine [%] | |", "|---|---|---|---|---|---|---|---|"]
    if all(v[2] is None for v in rows.values()):
        # Two levels only (fine grid not run): report the coarse -> medium change
        lines = [f"### {title}\n", "| Quantity | coarse | medium | change coarse -> medium [%] |", "|---|---|---|---|"]
        for name, vals in rows.items():
            c, m = vals[0], vals[1]
            lines.append(f"| {name} | {c:.4g} | {m:.4g} | {100 * (m - c) / abs(m):+.2f} |")
        return "\n".join(lines) + "\n"
    for name, vals in rows.items():
        if any(v is None or np.isnan(v) for v in vals):
            lines.append(f"| {name} | " + " | ".join("-" if v is None else f"{v:.4g}" for v in vals) + " | | | | incomplete |")
            continue
        g = gci(vals, r)
        lines.append(
            f"| {name} | {g['coarse']:.4g} | {g['medium']:.4g} | {g['fine']:.4g} | "
            + (f"{g['p']:.2f}" if np.isfinite(g["p"]) else "-")
            + " | "
            + (f"{g['ext']:.4g}" if np.isfinite(g["ext"]) else "-")
            + f" | {g['gci_pct']:.2f} | {g['note']} |"
        )
    return "\n".join(lines) + "\n"


out = []
for case, fn, metrics in [("fin", "lines.json", fin_metrics), ("flare", "lines.json", flare_metrics)]:
    for model in ["SA", "SA-Edwards"]:
        ds = [load(os.path.join(run_dir("adflow", case, lev, model), fn)) for lev in ["coarse", "medium", "fine"]]
        if not any(ds):
            continue
        ms = [metrics(d) if d else None for d in ds]
        keys = next(m for m in ms if m).keys()
        rows = {k: [m[k] if m else None for m in ms] for k in keys}
        cells = {"fin": "coarse 0.34 M / medium 0.97 M cells", "flare": "coarse 0.22 M / medium 0.63 M cells"}[case]
        out.append(table(f"{case} - {model} (r = sqrt(2); {cells})", rows, R3D))

for model in ["SA", "SA-Edwards"]:
    dirs = {lev: run_dir("adflow", "plate", lev, model) for lev in ["coarse", "medium", "fine"]}
    ds = {lev: load(os.path.join(d, "bl.json")) for lev, d in dirs.items()}
    ms = [plate_metrics(ds[lev], dirs[lev]) if ds[lev] else None for lev in ["coarse", "medium", "fine"]]
    if not any(ms):
        continue
    keys = next(m for m in ms if m).keys()
    rows = {k: [m[k] if m else None for m in ms] for k in keys}
    lab = "coarse / medium / fine (7k / 27k / 109k cells)"
    out.append(table(f"2-D plate, x = -11 cm - {model} (r = 2; {lab})", rows, R2D))

txt = "\n".join(out)
print(txt)
open(os.path.join(RESULTS, "grid_study.md"), "w").write(txt)
