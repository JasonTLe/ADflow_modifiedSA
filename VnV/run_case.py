"""Run a 3-D validation case of Edwards & Chandra (AIAA 94-2275).

Usage (through run.sh):  ./run.sh NPROC run_case.py {fin,flare} MODEL [restartFile]

  fin    Kussoy & Horstman Mach 8.2 flat plate / 15 deg sharp fin
  flare  Wideman et al. Mach 2.89 cylinder / offset 20 deg flare

MODEL is a key of common.MODELS; the grid is set with GRID=<file in grids/>.
Results are written to results/adflow/<case>/<level>/<MODEL>/ (level from the
grid name: medium, coarse or fine; or results/adflow/<case>/sensitivity/
<MODEL><RUN_TAG> if RUN_TAG is set): volume/surface solutions and the
extracted line data (lines.json)."""

import os
import sys
import glob
import json
import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from common import HERE, grid_level, run_dir, kussoy_ap, wideman_ap, solver_options, read_tecplot_surface, face_centers, gather_cells, winf, boundary_layer, init_boundary_layer, solve, dump_cells, init_from_cells

case = sys.argv[1]
model = sys.argv[2]
restart = sys.argv[3] if len(sys.argv) > 3 else None
comm = MPI.COMM_WORLD

if case == "fin":
    gridName = os.environ.get("GRID", "sharp_fin_M82.cgns")
    ap = kussoy_ap("fin")
elif case == "flare":
    gridName = os.environ.get("GRID", "cylinder_flare_M289.cgns")
    ap = wideman_ap("flare")
else:
    raise ValueError(case)
grid = os.path.join(HERE, "grids", gridName)
tag = os.environ.get("RUN_TAG", "")
out = run_dir("adflow", case, "sensitivity", model + tag) if tag else run_dir("adflow", case, grid_level(gridName), model)
os.makedirs(out, exist_ok=True)

extra = {"restartFile": restart} if restart else {}
opts = solver_options(grid, out, model, **extra)
CFDSolver = ADFLOW(options=opts, comm=comm)

interp = os.environ.get("INIT_FROM")  # .npz written by dump_cells (e.g. coarser grid)
if interp:
    CFDSolver.setAeroProblem(ap)
    init_from_cells(CFDSolver, ap, interp)
elif restart is None and os.environ.get("INIT_BL", "1") == "1":
    # Start from an approximate turbulent boundary layer instead of freestream
    if case == "fin":
        xLE, tanA = -1.64, np.tan(np.radians(15.0))

        def walls(xc):
            sp = np.maximum(xc[:, 0] - xLE, 0.0)
            sf = np.maximum(xc[:, 0], 0.0)
            dfin = np.where(xc[:, 0] > 0, (xc[:, 2] - xc[:, 0] * tanA) * np.cos(np.radians(15.0)), np.inf)
            return [(xc[:, 1], 0.028 * (sp / 1.53) ** 0.8), (dfin, 0.028 * (sf / 1.53) ** 0.8)]

        Tw = 300.0
    else:
        from meshes import flare_radius

        xLE = -0.90

        def walls(xc):
            r = np.hypot(xc[:, 1], xc[:, 2])
            phi = np.arctan2(xc[:, 2], xc[:, 1])
            d = r - flare_radius(xc[:, 0], phi)
            sp = np.maximum(xc[:, 0] - xLE, 0.0)
            return [(d, 0.011 * (sp / 0.88) ** 0.8)]

        Tw = None
    CFDSolver.setAeroProblem(ap)
    init_boundary_layer(CFDSolver, ap, walls, Twall=Tw)

solve(CFDSolver, ap, default_rk=1500 if (case == "fin" and restart is None and not interp) else 0)
dump_cells(CFDSolver, os.path.join(out, "cells.npz"))

res = {"case": case, "model": model, "converged": not ap.solveFailed}
comm.Barrier()

# Incoming boundary layer: fin case at x = -11 cm (far from the fin),
# flare case at x = -2 cm (phi = 90 deg, radial profile).
xc, w = gather_cells(CFDSolver)
wi = winf(CFDSolver)
if comm.rank == 0:
    if case == "fin":
        xs = -0.11
        xu = np.unique(np.round(xc[:, 0], 7))
        xi = xu[np.argmin(np.abs(xu - xs))]
        m = np.abs(xc[:, 0] - xi) < 1e-6
        zu = np.unique(np.round(xc[m, 2], 7))
        zi = zu[len(zu) // 2]
        m &= np.abs(xc[:, 2] - zi) < 1e-6
        yw = xc[m, 1]
    else:
        xs = -0.02
        xu = np.unique(np.round(xc[:, 0], 7))
        xi = xu[np.argmin(np.abs(xu - xs))]
        phi = np.degrees(np.arctan2(xc[:, 2], xc[:, 1]))
        m = np.abs(xc[:, 0] - xi) < 1e-6
        pu = np.unique(np.round(phi[m], 4))
        pi_ = pu[np.argmin(np.abs(pu - 90.0))]
        m &= np.abs(phi - pi_) < 1e-3
        yw = np.hypot(xc[m, 1], xc[m, 2]) - 0.0254
    d99, ds, th = boundary_layer(yw, w[m, 0], w[m, 1], wi[0], wi[1])
    res["incoming_BL"] = {"x": float(xi), "delta99_cm": d99 * 100, "deltaStar_cm": ds * 100, "theta_cm": th * 100}


def p_over_pinf(cp):
    return 1.0 + 0.5 * 1.4 * ap.mach**2 * cp


if comm.rank == 0:
    surf = read_tecplot_surface(sorted(glob.glob(os.path.join(out, "*_surf.plt")), key=os.path.getmtime)[-1])

    def lines(zone, var, axis, value, coord):
        """Element-centred values of var on the grid line whose element
        centres are closest to <axis> = value, sorted by <coord>."""
        X, v = face_centers(zone, var)
        a = X[:, axis]
        u = np.unique(np.round(a, 7))
        m = np.abs(a - u[np.argmin(np.abs(u - value))]) < 1e-6
        return X[m], v[m], float(u[np.argmin(np.abs(u - value))])

    if case == "fin":
        tanA = np.tan(np.radians(15.0))
        plate, fin = surf["plate"], surf["fin"]
        # Reference heat flux: undisturbed plate at the survey station x = -11 cm
        X, st = face_centers(plate, "StantonNumber")
        i = np.argmin(np.abs(X[:, 0] + 0.11) + np.abs(X[:, 2] - X[:, 2].max()))
        st_ref = st[i]
        res["St_ref"] = float(st_ref)
        for name, var, xs in [("plate_p_X18.2", "CoefPressure", 0.182), ("plate_q_X16.45", "StantonNumber", 0.1645), ("plate_cf_X15.5", "SkinFrictionMagnitude", 0.155)]:
            Xl, v, xa = lines(plate, var, 0, xs, 2)
            Z = (Xl[:, 2] - xa * tanA) * 100  # spanwise distance from the fin [cm]
            o = np.argsort(Z)
            if var == "CoefPressure":
                v = p_over_pinf(v)
            elif var == "StantonNumber":
                v = v / st_ref
            res[name] = {"x_actual": xa, "Z_cm": Z[o].tolist(), "value": v[o].tolist()}
        for name, var, xs in [("fin_p_X17.72", "CoefPressure", 0.1772), ("fin_q_X16.13", "StantonNumber", 0.1613)]:
            Xl, v, xa = lines(fin, var, 0, xs, 1)
            Y = Xl[:, 1] * 100
            o = np.argsort(Y)
            v = p_over_pinf(v) if var == "CoefPressure" else v / st_ref
            res[name] = {"x_actual": xa, "Y_cm": Y[o].tolist(), "value": v[o].tolist()}
    else:
        body = surf["body"]
        X, cfx = face_centers(body, "SkinFrictionX")
        _, cp = face_centers(body, "CoefPressure")
        phi = np.degrees(np.arctan2(X[:, 2], X[:, 1]))
        for ph in [0, 45, 90, 135, 180]:
            u = np.unique(np.round(phi, 4))
            p0 = u[np.argmin(np.abs(u - ph))]
            m = np.abs(phi - p0) < 1e-3
            o = np.argsort(X[m, 0])
            xl = X[m, 0][o] * 100
            c = cfx[m][o]
            # Separation / attachment: sign changes of Cfx (upstream of the afterbody)
            s = np.where(np.diff(np.sign(c)) != 0)[0]
            crossings = [float(xl[k] - c[k] * (xl[k + 1] - xl[k]) / (c[k + 1] - c[k])) for k in s if -5 < xl[k] < 15]
            res[f"phi{ph}"] = {
                "phi_actual": float(p0),
                "x_cm": xl.tolist(),
                "Cfx": c.tolist(),
                "P_over_Pinf": p_over_pinf(cp[m][o]).tolist(),
                "Cfx_zero_crossings_cm": crossings,
                "Cfx_upstream_x-3cm": float(np.interp(-3.0, xl, c)),
            }
    json.dump(res, open(os.path.join(out, "lines.json"), "w"))
    summary = {k: v for k, v in res.items() if not isinstance(v, dict)}
    summary["incoming_BL"] = res["incoming_BL"]
    print("RESULT", json.dumps(summary))
