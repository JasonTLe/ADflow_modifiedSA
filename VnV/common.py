"""Shared setup for the Edwards & Chandra (AIAA 94-2275) validation cases."""

import os
import json
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem

HERE = os.path.dirname(os.path.abspath(__file__))

# Output layout: results/<code>/<case>/<level>/<model>, where <level> is the
# grid level (medium = the grid used for all reported results; coarse and fine
# only for the grid-refinement study) or "sensitivity" for the sensitivity
# checks. Documents (tables) go to results/, plots to results/figures/.
RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(RESULTS, "figures")
GRID_LEVELS = {"": "medium", "_refc": "coarse", "_reff": "fine"}


def grid_level(gridName):
    """Grid level from a grid file name (..._refc.cgns = coarse, ..._reff.cgns = fine)."""
    for tag, level in GRID_LEVELS.items():
        if tag and gridName.replace(".cgns", "").endswith(tag):
            return level
    return "medium"


def run_dir(code, case, level, model):
    """results/<code>/<case>/<level>/<model>"""
    return os.path.join(RESULTS, code, case, level, model)

# Turbulence model variants to compare. Keys are used in file names.
MODELS = {
    "SA": {"turbulenceModel": "SA", "useft2SA": False},
    "SA-Edwards": {"turbulenceModel": "SA-Edwards"},
    "SA-comp": {"turbulenceModel": "SA", "useft2SA": False, "useCompressibilitySA": True},
}


def kussoy_ap(name="kussoy"):
    """Kussoy & Horstman 15 deg sharp fin, Table in Sec. 6.1 of Edwards & Chandra."""
    return AeroProblem(
        name=name,
        mach=8.20,
        T=76.89,
        reynolds=5.32e6,
        reynoldsLength=1.0,
        alpha=0.0,
        beta=0.0,
        areaRef=1.0,
        chordRef=1.0,
        R=287.05,
        evalFuncs=["cd"],
    )


def wideman_ap(name="wideman"):
    """Wideman et al. cylinder / offset flare, Table 1 of NASA CR-202617
    (P = 5.54 kPa, T = 105 K, M = 2.89; Re = 15e6 / m)."""
    return AeroProblem(
        name=name,
        mach=2.89,
        T=105.0,
        P=5540.0,
        alpha=0.0,
        beta=0.0,
        areaRef=1.0,
        chordRef=1.0,
        R=287.05,
        evalFuncs=["cd"],
    )


def solver_options(gridFile, outDir, model, **extra):
    opts = {
        "gridFile": gridFile,
        "outputDirectory": outDir,
        "equationType": "RANS",
        # Upwind (Roe) scheme for the high Mach number flow
        "discretization": "upwind",
        "limiter": "van Albada",
        "smoother": "DADI",
        "MGCycle": "sg",
        "useApproxWallDistance": False,
        "CFL": 1.0,
        "nCycles": 20000,
        "nSubiterTurb": 5,
        "L2Convergence": 1e-9,
        "useANKSolver": True,
        # At Mach 8 the ANK solver cannot get through the start-up
        # transient with these settings; there an explicit Runge-Kutta
        # phase is run first (STAGE1_RK, see solve()).
        "ANKSwitchTol": 1e3,
        "ANKCFLLimit": 1e4,
        "ANKSecondOrdSwitchTol": 1e-3,
        "ANKCoupledSwitchTol": 1e-5,
        "useNKSolver": False,
        "monitorVariables": ["resrho", "resturb", "cd"],
        "surfaceVariables": ["cp", "p", "temp", "cf", "cfx", "cfy", "cfz", "ch", "yplus", "vx", "vy", "vz"],
        "volumeVariables": ["resrho", "resturb", "mach", "eddy", "temp"],
        "writeTecplotSurfaceSolution": True,
        "writeVolumeSolution": True,
        "printIterations": True,
        "printAllOptions": False,
    }
    opts.update(MODELS[model])
    opts.update(extra)
    # Optional overrides for solver-setting studies, e.g. ADFLOW_EXTRA='{"ANKCFL0": 1.0}'
    opts.update(json.loads(os.environ.get("ADFLOW_EXTRA", "{}")))
    return opts


# ------------------------------------------------------------------------------
# Post-processing
# ------------------------------------------------------------------------------
def gather_cells(CFDSolver):
    """Return (xCen[n,3], w[n,6]) of all cells on rank 0 (None elsewhere).
    w are ADflow's nondimensional states (rho, u, v, w, rhoE, nuTilde)."""
    comm = CFDSolver.comm
    nw = 6
    w = CFDSolver.getStates().reshape(-1, nw)
    x = CFDSolver.adflow.utils.getcellcenters(1, w.shape[0]).T
    xs = comm.gather(x, root=0)
    ws = comm.gather(w, root=0)
    if comm.rank == 0:
        return np.vstack(xs), np.vstack(ws)
    return None, None


def winf(CFDSolver):
    return np.array(CFDSolver.adflow.flowvarrefstate.winf[:6])


def boundary_layer(y, rho, u, rhoInf, uInf):
    """delta99, displacement and momentum thickness of a profile u(y)."""
    o = np.argsort(y)
    y, rho, u = y[o], rho[o], u[o]
    ue = uInf
    i99 = np.argmax(u >= 0.99 * ue)
    d99 = np.interp(0.99 * ue, u[max(i99 - 1, 0) : i99 + 1], y[max(i99 - 1, 0) : i99 + 1]) if i99 > 0 else np.nan
    f1 = 1 - rho * u / (rhoInf * uInf)
    f2 = rho * u / (rhoInf * uInf) * (1 - u / uInf)
    mask = y <= 1.5 * d99
    yy = np.concatenate([[0.0], y[mask]])
    ds = np.trapz(np.concatenate([[1.0], f1[mask]]), yy)
    th = np.trapz(np.concatenate([[0.0], f2[mask]]), yy)
    return d99, ds, th


def read_tecplot_surface(fileName):
    """Reader for ADflow's binary Tecplot surface files (#!TDV112), following
    writeTecplotSurfaceFile in src/output/tecplotIO.F90. Returns a dict
    {familyName: {"nodes": dict of node arrays, "conn": (nElem, 4) array}}.
    Zones share the node data (all surface nodes) and differ by connectivity."""
    buf = open(fileName, "rb").read()
    pos = [8]  # skip '#!TDV112'

    def i32(n=1):
        v = np.frombuffer(buf, "<i4", n, pos[0])
        pos[0] += 4 * n
        return v if n > 1 else int(v[0])

    def f32(n=1):
        v = np.frombuffer(buf, "<f4", n, pos[0])
        pos[0] += 4 * n
        return v if n > 1 else float(v[0])

    def f64(n=1):
        v = np.frombuffer(buf, "<f8", n, pos[0])
        pos[0] += 8 * n
        return v if n > 1 else float(v[0])

    def string():
        chars = []
        while True:
            c = i32()
            if c == 0:
                return "".join(chars)
            chars.append(chr(c))

    i32(), i32()  # byte order, file type
    string()  # title
    nv = i32()
    names = [string() for _ in range(nv)]
    zones = []
    while True:
        marker = f32()
        if abs(marker - 357.0) < 1e-3:
            break
        name = string()
        i32(), i32()  # parent, strand
        f64()  # solution time
        i32()  # color
        ztype, packing, varloc, faceNb = i32(), i32(), i32(), i32()
        nNodes, nElem = i32(), i32()
        i32(), i32(), i32(), i32()  # cell dims, aux
        zones.append({"name": name, "nNodes": nNodes, "nElem": nElem})
    shared = None
    for z in zones:
        f32()  # zone marker
        dtypes = i32(nv)
        i32()  # passive
        share = i32()
        if share == 0:
            i32()  # connectivity sharing
            f64(2 * nv)  # min/max
            data = {}
            for iv, nm in enumerate(names):
                data[nm] = f32(z["nNodes"]) if dtypes[iv] == 1 else f64(z["nNodes"])
            shared = data
        else:
            i32(nv)
            i32()
            data = shared
        conn = i32(4 * z["nElem"]).reshape(-1, 4)
        z["nodes"] = data
        z["conn"] = conn
    return {z["name"]: z for z in zones}


def face_centers(zone, var):
    """Element-centred coordinates and element-averaged variable."""
    c = zone["conn"]
    n = zone["nodes"]
    X = np.stack([n["CoordinateX"], n["CoordinateY"], n["CoordinateZ"]], 1)
    return X[c].mean(1), n[var][c].mean(1)


# ------------------------------------------------------------------------------
# Initial condition: approximate compressible turbulent boundary layer
# ------------------------------------------------------------------------------
def init_boundary_layer(CFDSolver, ap, walls, Twall=None, recovery=0.89):
    """Overwrite the (freestream) initial states with an approximate turbulent
    boundary layer to avoid the impulsive-start transient at high Mach number.

    walls: function xc[n,3] -> list of (d, delta) arrays, one pair per wall:
           wall distance d and local boundary-layer thickness delta (m).
    The velocity is (y/delta)^(1/7) (product over walls) along the freestream
    direction, the temperature follows the Crocco-Busemann relation (to Twall,
    or to the adiabatic wall temperature if Twall is None) at constant
    pressure, and nuTilde follows a mixing-length estimate."""
    fv = CFDSolver.adflow.flowvarrefstate
    gam = float(fv.gammainf)
    w0 = CFDSolver.getStates().reshape(-1, 6).copy()
    n = w0.shape[0]
    xc = CFDSolver.adflow.utils.getcellcenters(1, n).T
    wi = np.array(fv.winf[:6])
    pinf = float(fv.pinf)
    rhoinf = wi[0]
    V = np.linalg.norm(wi[1:4])
    M = ap.mach
    Tinf = pinf / rhoinf  # nondimensional (R = 1 in these units up to a constant)
    Taw = Tinf * (1 + recovery * 0.5 * (gam - 1) * M**2)
    Tw = Taw if Twall is None else Twall / ap.T * Tinf

    f = np.ones(n)
    mixing = np.full(n, np.inf)
    for d, delta in walls(xc):
        eta = np.clip(d / np.maximum(delta, 1e-12), 0.0, 1.0)
        f *= eta ** (1.0 / 7.0)
        # nu_t ~ kappa u_tau y (1 - y/delta), u_tau ~ 0.04 U
        mixing = np.minimum(mixing, 0.41 * 0.04 * V * d * (1.0 - eta) + np.where(eta >= 1, np.inf, 0.0))
    uu = f  # u / U
    T = Tw + (Taw - Tw) * uu - (Taw - Tinf) * uu**2
    rho = pinf / T
    w = w0.copy()
    w[:, 0] = rho
    w[:, 1:4] = wi[1:4][None, :] * uu[:, None]
    w[:, 4] = pinf / (gam - 1) + 0.5 * rho * (V * uu) ** 2
    nut = np.where(np.isfinite(mixing), mixing, 0.0) * rho / rhoinf
    w[:, 5] = np.maximum(wi[5], nut * rhoinf / rho) if True else wi[5]
    CFDSolver.setStates(w.flatten())
    return w


def solve(CFDSolver, ap, default_rk=0):
    """Solve, optionally with an explicit Runge-Kutta start-up phase
    (STAGE1_RK=<n> iterations at CFL 1) before the ANK solver takes over
    with its standard CFL ramp."""
    nRK = int(os.environ.get("STAGE1_RK", str(default_rk)))
    if nRK > 0:
        saved = {k: CFDSolver.getOption(k) for k in ["useANKSolver", "smoother", "CFL", "nCycles"]}
        CFDSolver.setOption("useANKSolver", False)
        CFDSolver.setOption("smoother", "Runge-Kutta")
        CFDSolver.setOption("CFL", 1.0)
        CFDSolver.setOption("nCycles", nRK)
        CFDSolver(ap, writeSolution=False)
        for k, v in saved.items():
            CFDSolver.setOption(k, v)
        for k, v in json.loads(os.environ.get("STAGE2_EXTRA", "{}")).items():
            CFDSolver.setOption(k, v)
    CFDSolver(ap)


# ------------------------------------------------------------------------------
# Initial condition interpolated from a solution on another grid
# ------------------------------------------------------------------------------
def dump_cells(CFDSolver, fileName):
    """Save all cell centres and states (rank 0) to an .npz file."""
    x, w = gather_cells(CFDSolver)
    if CFDSolver.comm.rank == 0:
        np.savez(fileName, x=x, w=w)


def init_from_cells(CFDSolver, ap, fileName, k=8):
    """Set the states by inverse-distance interpolation (k nearest cell
    centres) of a solution saved with dump_cells, e.g. from a coarser grid of
    the same family. Used only as an initial condition."""
    from scipy.spatial import cKDTree

    d = np.load(fileName)
    tree = cKDTree(d["x"])
    n = CFDSolver.getStateSize() // 6
    xc = CFDSolver.adflow.utils.getcellcenters(1, n).T
    dist, idx = tree.query(xc, k=k)
    wgt = 1.0 / np.maximum(dist, 1e-14) ** 2
    wgt /= wgt.sum(1, keepdims=True)
    w = np.einsum("nk,nkv->nv", wgt, d["w"][idx])
    CFDSolver.setStates(w.flatten())
