"""Extract line data from the ADflow surface solution of a validation case.

Usage: python postprocess.py {fin,flare} LEVEL MODEL   (LEVEL: medium, coarse, fine)
Reads results/adflow/<case>/<LEVEL>/<MODEL>/*_surf.plt and writes lines.json there
(the incoming boundary layer data written by run_case.py is kept)."""

import os
import sys
import glob
import json
import numpy as np
from common import HERE, read_tecplot_surface, run_dir

MACH = {"fin": 8.2, "flare": 2.89}


def nodes(zone):
    n = zone["nodes"]
    X = np.stack([n["CoordinateX"], n["CoordinateY"], n["CoordinateZ"]], 1)
    # Surface nodes are duplicated across blocks; keep unique coordinates
    _, idx = np.unique(np.round(X, 9), axis=0, return_index=True)
    return X[idx], {k: v[idx] for k, v in n.items()}


def p_over_pinf(cp, mach):
    return 1.0 + 0.5 * 1.4 * mach**2 * cp


def grid_line(a, value, tol=1e-7):
    """Mask of nodes on the grid line a = const closest to value."""
    u = np.unique(np.round(a, 7))
    a0 = u[np.argmin(np.abs(u - value))]
    return np.abs(a - a0) < tol + 1e-9, float(a0)


def post_fin(surf):
    tanA = np.tan(np.radians(15.0))
    res = {}
    Xp, vp = nodes(surf["plate"])
    Xf, vf = nodes(surf["fin"])
    # Reference heat flux: undisturbed plate at the survey station x = -11 cm, far from the fin
    m, xa = grid_line(Xp[:, 0], -0.11)
    st_ref = float(np.median(vp["StantonNumber"][m]))
    res["St_ref"] = st_ref
    for name, var, xs in [
        ("plate_p_X18.2", "CoefPressure", 0.182),
        ("plate_q_X16.45", "StantonNumber", 0.1645),
        ("plate_cf_X15.5", "SkinFrictionMagnitude", 0.155),
    ]:
        m, xa = grid_line(Xp[:, 0], xs)
        Z = (Xp[m, 2] - xa * tanA) * 100  # spanwise distance from the fin [cm]
        v = vp[var][m]
        # Drop the plate/fin corner-line node (nodal average of both walls)
        keep = Z > 1e-4
        Z, v = Z[keep], v[keep]
        v = p_over_pinf(v, 8.2) if var == "CoefPressure" else (v / st_ref if var == "StantonNumber" else v)
        o = np.argsort(Z)
        res[name] = {"x_actual": xa, "Z_cm": Z[o].tolist(), "value": v[o].tolist()}
    for name, var, xs in [("fin_p_X17.72", "CoefPressure", 0.1772), ("fin_q_X16.13", "StantonNumber", 0.1613)]:
        m, xa = grid_line(Xf[:, 0], xs)
        Y = Xf[m, 1] * 100
        v = vf[var][m]
        # Drop the corner-line node and the upper part influenced by the
        # fin tip at the domain top (y = 15 cm)
        keep = (Y > 1e-4) & (Y < 12.0)
        Y, v = Y[keep], v[keep]
        v = p_over_pinf(v, 8.2) if var == "CoefPressure" else v / st_ref
        o = np.argsort(Y)
        res[name] = {"x_actual": xa, "Y_cm": Y[o].tolist(), "value": v[o].tolist()}
    return res


def post_flare(surf):
    res = {}
    X, v = nodes(surf["body"])
    phi = np.degrees(np.arctan2(X[:, 2], X[:, 1]))
    for ph in [0, 45, 90, 135, 180]:
        m = np.abs(phi - ph) < 0.05 if ph < 180 else np.abs(np.abs(phi) - 180) < 0.05
        o = np.argsort(X[m, 0])
        xl = X[m, 0][o] * 100
        c = v["SkinFrictionX"][m][o]
        P = p_over_pinf(v["CoefPressure"][m][o], 2.89)
        s = np.where(np.diff(np.sign(c)) != 0)[0]
        crossings = [float(xl[k] - c[k] * (xl[k + 1] - xl[k]) / (c[k + 1] - c[k])) for k in s if -10 < xl[k] < 15]
        # Undisturbed values well upstream of the interaction
        up = (xl > -12) & (xl < -8)
        res[f"phi{ph}"] = {
            "x_cm": xl.tolist(),
            "Cfx": c.tolist(),
            "P_over_Pinf": P.tolist(),
            "Cfx_zero_crossings_cm": crossings,
            "Cfx_upstream_x-10cm": float(np.mean(c[up])) if up.any() else None,
            "P_upstream_x-10cm": float(np.mean(P[up])) if up.any() else None,
        }
    return res


if __name__ == "__main__":
    case, level, model = sys.argv[1], sys.argv[2], sys.argv[3]
    out = run_dir("adflow", case, level, model)
    surf = read_tecplot_surface(sorted(glob.glob(os.path.join(out, "*_surf.plt")), key=os.path.getmtime)[-1])
    res = post_fin(surf) if case == "fin" else post_flare(surf)
    fn = os.path.join(out, "lines.json")
    old = json.load(open(fn)) if os.path.exists(fn) else {}
    for k in ["case", "model", "converged", "incoming_BL"]:
        if k in old:
            res[k] = old[k]
    json.dump(res, open(fn, "w"))
    summary = {k: v for k, v in res.items() if not isinstance(v, (dict, list))}
    if case == "flare":
        for ph in [0, 90, 180]:
            r = res[f"phi{ph}"]
            summary[f"phi{ph}"] = {"Cfx_up": r["Cfx_upstream_x-10cm"], "P_up": r["P_upstream_x-10cm"], "zero_crossings": [round(c, 2) for c in r["Cfx_zero_crossings_cm"]]}
    print(json.dumps(summary, indent=1))
