"""SU2 version of the 2-D Mach 8.2 cold-wall flat plate, on exactly the same
grids as ADflow (VnV/grids/flat_plate_M82_LE164*.cgns, the z = 0 plane).

Usage: python make_su2_case.py LEVEL MODEL
  LEVEL: coarse | medium | fine
  MODEL: SA | SA-Edwards
Writes results/su2/plate/<LEVEL>/<MODEL>/{plate.su2, stage1.cfg, stage2.cfg}.

The physical modelling is matched to ADflow's defaults: SA without ft2,
Sutherland (mu_ref 1.716e-5, T_ref 273.15 K, S 110.55 K), Pr 0.72,
Pr_t 0.9, gamma 1.4, R 287.05, freestream nuTilde = 1.342 nu (ADflow's
eddyVisInfRatio 0.009), isothermal wall 300 K, Roe + MUSCL/van Albada for
the flow, first-order upwind for the turbulence."""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from cgnsutilities.cgnsutilities import readGrid  # noqa: E402
from common import run_dir, GRID_LEVELS  # noqa: E402

level = sys.argv[1] if len(sys.argv) > 1 else "medium"
model = sys.argv[2] if len(sys.argv) > 2 else "SA"
tag = {v: k for k, v in GRID_LEVELS.items()}[level]
out = run_dir("su2", "plate", level, model)
os.makedirs(out, exist_ok=True)

g = readGrid(os.path.join(HERE, "..", "grids", f"flat_plate_M82_LE164{tag}.cgns"))
B = {b.name: b for b in g.blocks}
lead = B["leadin"].coords[:, :, 0, :2]
plate = B["plate"].coords[:, :, 0, :2]
assert np.allclose(lead[-1], plate[0])
X = np.concatenate([lead, plate[1:]], axis=0)  # (ni, nj, 2)
ni, nj = X.shape[:2]
nLead = lead.shape[0]
idx = np.arange(ni * nj).reshape(ni, nj)  # node id = i*nj + j

with open(os.path.join(out, "plate.su2"), "w") as f:
    f.write("NDIME= 2\n")
    f.write(f"NELEM= {(ni - 1) * (nj - 1)}\n")
    e = 0
    for i in range(ni - 1):
        for j in range(nj - 1):
            f.write(f"9 {idx[i, j]} {idx[i + 1, j]} {idx[i + 1, j + 1]} {idx[i, j + 1]} {e}\n")
            e += 1
    f.write(f"NPOIN= {ni * nj}\n")
    for i in range(ni):
        for j in range(nj):
            f.write(f"{X[i, j, 0]:.15e} {X[i, j, 1]:.15e} {idx[i, j]}\n")
    markers = {
        "inlet": [(idx[0, j], idx[0, j + 1]) for j in range(nj - 1)],
        "outlet": [(idx[-1, j], idx[-1, j + 1]) for j in range(nj - 1)],
        "top": [(idx[i, -1], idx[i + 1, -1]) for i in range(ni - 1)],
        "symm": [(idx[i, 0], idx[i + 1, 0]) for i in range(nLead - 1)],
        "wall": [(idx[i, 0], idx[i + 1, 0]) for i in range(nLead - 1, ni - 1)],
    }
    f.write(f"NMARK= {len(markers)}\n")
    for name, edges in markers.items():
        f.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(edges)}\n")
        for a, b in edges:
            f.write(f"3 {a} {b}\n")

sa_opts = {"SA": "NONE", "SA-Edwards": "EDWARDS"}[model]
CFG = """% SU2 Mach 8.2 cold-wall flat plate (Kussoy & Horstman inflow), model {model}
SOLVER= RANS
KIND_TURB_MODEL= SA
SA_OPTIONS= {sa_opts}
MATH_PROBLEM= DIRECT
RESTART_SOL= {restart}
READ_BINARY_RESTART= NO

MACH_NUMBER= 8.2
AOA= 0.0
SIDESLIP_ANGLE= 0.0
INIT_OPTION= REYNOLDS
FREESTREAM_OPTION= TEMPERATURE_FS
FREESTREAM_TEMPERATURE= 76.89
REYNOLDS_NUMBER= 5.32E6
REYNOLDS_LENGTH= 1.0
FREESTREAM_NU_FACTOR= 1.342
REF_DIMENSIONALIZATION= DIMENSIONAL

FLUID_MODEL= STANDARD_AIR
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.05
VISCOSITY_MODEL= SUTHERLAND
MU_REF= 1.716E-5
MU_T_REF= 273.15
SUTHERLAND_CONSTANT= 110.55
PRANDTL_LAM= 0.72
PRANDTL_TURB= 0.90

REF_ORIGIN_MOMENT_X= 0.0
REF_ORIGIN_MOMENT_Y= 0.0
REF_ORIGIN_MOMENT_Z= 0.0
REF_LENGTH= 1.0
REF_AREA= 1.0

MARKER_ISOTHERMAL= ( wall, 300.0 )
MARKER_SYM= ( symm )
MARKER_FAR= ( inlet, top )
MARKER_SUPERSONIC_OUTLET= ( outlet )
MARKER_PLOTTING= ( wall )
MARKER_MONITORING= ( wall )

NUM_METHOD_GRAD= GREEN_GAUSS
CONV_NUM_METHOD_FLOW= ROE
ENTROPY_FIX_COEFF= 0.0
MUSCL_FLOW= {muscl}
SLOPE_LIMITER_FLOW= VAN_ALBADA_EDGE
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_FLOW= EULER_IMPLICIT
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_NUMBER= {cfl}
CFL_ADAPT= NO
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 20

ITER= {niter}
CONV_FIELD= RMS_DENSITY
CONV_RESIDUAL_MINVAL= -11

MESH_FILENAME= plate.su2
MESH_FORMAT= SU2
OUTPUT_FILES= ( RESTART_ASCII, PARAVIEW_ASCII )
VOLUME_FILENAME= flow
SOLUTION_FILENAME= restart_flow.csv
RESTART_FILENAME= restart_flow.csv
SURFACE_FILENAME= surface_flow
CONV_FILENAME= history
SCREEN_OUTPUT= ( INNER_ITER, RMS_DENSITY, RMS_NU_TILDE, CUR_TIME )
HISTORY_OUTPUT= ( ITER, RMS_RES )
VOLUME_OUTPUT= ( COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL )
OUTPUT_WRT_FREQ= {wrt}
"""
# Stage 1: first order, stage 2: second order (MUSCL + van Albada) restarted
# from stage 1. The corner node at the plate leading edge keeps SU2's global
# residual from dropping (see README); convergence is judged by the residual
# away from the leading edge and by the stationarity of the wall quantities.
open(os.path.join(out, "stage1.cfg"), "w").write(CFG.format(model=model, sa_opts=sa_opts, cfl=5.0, muscl="NO", restart="NO", niter=8000, wrt=8000))
open(os.path.join(out, "stage2.cfg"), "w").write(CFG.format(model=model, sa_opts=sa_opts, cfl=2.0, muscl="YES", restart="YES", niter=12000, wrt=2000) + "WRT_VOLUME_OVERWRITE= NO\n")

print(out, ni, nj)
