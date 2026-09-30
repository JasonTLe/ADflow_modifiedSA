"""2-D Mach 8.2 cold-wall flat plate: incoming boundary layer of the Kussoy &
Horstman sharp-fin experiment. The fin apex is at x = 0, 176 cm downstream of
the plate leading edge; the boundary layer was surveyed at x = -11 cm
(delta = 3.7 cm, theta = 0.094 cm, Cf = 1.0e-3, qw = 1.04 W/cm^2)."""

import os
import sys
import json
import numpy as np
from mpi4py import MPI
from adflow import ADFLOW
from common import HERE, grid_level, run_dir, kussoy_ap, solver_options, gather_cells, winf, boundary_layer, init_boundary_layer, solve

model = sys.argv[1] if len(sys.argv) > 1 else "SA"
gridName = sys.argv[2] if len(sys.argv) > 2 else "flat_plate_M82_LE164.cgns"
# results/adflow/plate/<level>/<model>, or .../sensitivity/<model><RUN_TAG>
tag = os.environ.get("RUN_TAG", "")
out = run_dir("adflow", "plate", "sensitivity", model + tag) if tag else run_dir("adflow", "plate", grid_level(gridName), model)
os.makedirs(out, exist_ok=True)

opts = solver_options(os.path.join(HERE, "grids", gridName), out, model)
ap = kussoy_ap("flatplate")

CFDSolver = ADFLOW(options=opts, comm=MPI.COMM_WORLD)
if os.environ.get("INIT_BL", "1") == "1":
    xLE = float(os.environ.get("XLE", "-1.64"))

    def walls(xc):
        s = np.maximum(xc[:, 0] - xLE, 0.0)
        return [(xc[:, 1], 0.028 * (s / 1.53) ** 0.8)]

    CFDSolver.setAeroProblem(ap)
    init_boundary_layer(CFDSolver, ap, walls, Twall=300.0)
solve(CFDSolver, ap, default_rk=1500)

x, w = gather_cells(CFDSolver)
wi = winf(CFDSolver)
if MPI.COMM_WORLD.rank == 0:
    res = {"model": model, "converged": not ap.solveFailed}
    zc = np.unique(np.round(x[:, 2], 8))[0]
    for xs in [-0.11, -0.05, 0.0]:
        sel = np.abs(x[:, 2] - zc) < 1e-9
        xu = np.unique(np.round(x[sel, 0], 8))
        xi = xu[np.argmin(np.abs(xu - xs))]
        m = sel & (np.abs(x[:, 0] - xi) < 1e-7)
        d99, ds, th = boundary_layer(x[m, 1], w[m, 0], w[m, 1], wi[0], wi[1])
        res[f"x={xs:+.2f}"] = {"x": float(xi), "delta99_cm": d99 * 100, "deltaStar_cm": ds * 100, "theta_cm": th * 100}
    # Wall quantities from the surface file
    from common import read_tecplot_surface, face_centers
    import glob
    z = read_tecplot_surface(sorted(glob.glob(os.path.join(out, "*_surf.plt")), key=os.path.getmtime)[-1])["plate"]
    X, cf = face_centers(z, "SkinFrictionMagnitude")
    _, st = face_centers(z, "StantonNumber")
    T0 = 76.89 * (1 + 0.2 * 8.2**2)
    rhoU_cp_dT = ap.rho * ap.V * 1004.5 * (T0 - 300.0)  # W/m^2 per unit Stanton number
    for key in list(res):
        if key.startswith("x="):
            i = np.argmin(np.abs(X[:, 0] - res[key]["x"]))
            res[key]["Cf"] = float(cf[i])
            res[key]["qw_W_cm2"] = float(st[i] * rhoU_cp_dT / 1e4)
    print("RESULT", json.dumps(res))
    json.dump(res, open(os.path.join(out, "bl.json"), "w"), indent=2)
