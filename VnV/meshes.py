"""
Structured grids for the validation cases of Edwards & Chandra, AIAA 94-2275:

* flat_plate_2d:     flat plate (2 cells thick), used to check the incoming
                     boundary layer of the Kussoy-Horstman experiment.
* sharp_fin:         Kussoy & Horstman Mach 8.2 flat plate / 15 deg sharp fin.
* cylinder_flare:    Wideman et al. Mach 2.89 cylinder / offset 20 deg flare.

All lengths are in meters. Written with cgnsutilities.
"""

import numpy as np
from cgnsutilities.cgnsutilities import Grid, Block, BocoDataSet, BocoDataSetArray, CGNSDATATYPES


# ------------------------------------------------------------------------------
# 1-D point distributions
# ------------------------------------------------------------------------------
def tanh_one_sided(n, length, ds0):
    """n points on [0, length] with first spacing ds0 (hyperbolic tangent)."""
    from scipy.optimize import brentq

    xi = np.linspace(0.0, 1.0, n)

    def first(delta):
        s = 1.0 + np.tanh(delta * (xi - 1.0)) / np.tanh(delta)
        return s[1] * length - ds0

    if first(1e-6) < 0:  # uniform spacing already small enough
        return xi * length
    delta = brentq(first, 1e-6, 50.0)
    return (1.0 + np.tanh(delta * (xi - 1.0)) / np.tanh(delta)) * length


def geometric_to(n, length, ds0):
    """n points on [0, length], first spacing ds0, geometric growth."""
    from scipy.optimize import brentq

    if abs(ds0 * (n - 1) - length) < 1e-14:
        return np.linspace(0, length, n)
    f = lambda r: ds0 * (r ** (n - 1) - 1.0) / (r - 1.0) - length
    r = brentq(f, 1.0 + 1e-10, 3.0)
    return np.concatenate([[0.0], np.cumsum(ds0 * r ** np.arange(n - 1))])


def piecewise(segments):
    """Concatenate (x0, x1, n, ds_start, ds_end) segments using a two-sided
    tanh-like distribution (Vinokur) approximated by blending."""
    pts = []
    for x0, x1, n, dsa, dsb in segments:
        L = x1 - x0
        s = vinokur(n, L, dsa, dsb) + x0
        if pts:
            s = s[1:]
        pts.append(s)
    return np.concatenate(pts)


def vinokur(n, L, ds1, ds2):
    """Two-sided Vinokur stretching on [0, L] with end spacings ds1, ds2."""
    from scipy.optimize import brentq

    xi = np.linspace(0, 1, n)
    A = np.sqrt(ds2 / ds1)
    B = L / ((n - 1) * np.sqrt(ds1 * ds2))
    if abs(B - 1.0) < 1e-6:
        u = xi
    elif B > 1.0:
        dy = brentq(lambda d: np.sinh(d) / d - B, 1e-8, 100)
        u = 0.5 * (1 + np.tanh(dy * (xi - 0.5)) / np.tanh(dy / 2))
    else:
        dy = brentq(lambda d: np.sin(d) / d - B, 1e-8, np.pi - 1e-8)
        u = 0.5 * (1 + np.tan(dy * (xi - 0.5)) / np.tan(dy / 2))
    s = u / (A + (1 - A) * u)
    return s * L


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
def _isothermal_dataset(Twall):
    ds = BocoDataSet("BCDataSet", "bcwallviscousisothermal")
    arr = np.array([Twall], dtype=np.float64)
    dims = np.ones(3, dtype=np.int32, order="F")
    dims[0] = 1
    ds.addDirichletDataSet(BocoDataSetArray("Temperature", CGNSDATATYPES["RealDouble"], 1, dims, arr))
    return [ds]


def _write(blocks, bcs, fileName):
    """blocks: list of (name, coords[ni,nj,nk,3]); bcs: list of
    (blockIndex, face, bcType, family, dataSet or None)."""
    grid = Grid()
    for name, X in blocks:
        X = np.asfortranarray(X)
        grid.addBlock(Block(name, np.array(X.shape[:3]), X))
    grid.connect()
    for ib, face, bcType, fam, ds in bcs:
        grid.blocks[ib].overwriteBCs(face, bcType, fam, ds if ds is not None else [])
    grid.writeToCGNS(fileName)
    return grid


def _n(intervals, ref):
    """Number of points for a segment of `intervals` intervals refined by ref."""
    return int(round(intervals * ref)) + 1


# ------------------------------------------------------------------------------
# Kussoy-Horstman flat plate / sharp fin (Mach 8.2)
# ------------------------------------------------------------------------------
def fin_x_distribution(xLE=-1.76, xEnd=0.25, ref=1.0):
    """Streamwise nodes: plate leading edge at xLE, fin apex at x=0.
    ref scales the number of intervals (and divides the spacings)."""
    return piecewise(
        [
            (xLE, -0.12, _n(41, ref), 2.0e-4 / ref, 1.5e-2 / ref),
            (-0.12, 0.0, _n(41, ref), 1.5e-2 / ref, 1.5e-3 / ref),
            (0.0, xEnd, _n(81, ref), 1.5e-3 / ref, 3.0e-3 / ref),
        ]
    )


def flat_plate_2d(fileName, xLE=-1.76, xEnd=0.25, yTop=0.15, ny=81, dy0=2.0e-5, Twall=300.0, xUp=0.02, ref=1.0):
    """2-D flat plate (2 cells thick in z) with a short inviscid lead-in."""
    xin = np.linspace(xLE - xUp, xLE, _n(8, ref))
    xp = fin_x_distribution(xLE, xEnd, ref)
    ny = _n(ny - 1, ref)
    y = tanh_one_sided(ny, yTop, dy0 / ref)
    z = np.array([0.0, 0.005, 0.01])
    blocks = []
    for name, x in [("leadin", xin), ("plate", xp)]:
        X = np.zeros((len(x), ny, 3, 3))
        X[..., 0] = x[:, None, None]
        X[..., 1] = y[None, :, None]
        X[..., 2] = z[None, None, :]
        blocks.append((name, X))
    bcs = [
        (0, "ilow", "bcinflowssupersonic", "inflow", None),
        (0, "jlow", "bcsymmetryplane", "symm_plate", None),
        (1, "jlow", "bcwallviscousisothermal", "plate", _isothermal_dataset(Twall)),
        (1, "ihigh", "bcoutflowsupersonic", "outflow", None),
    ]
    for ib in range(2):
        bcs += [
            (ib, "jhigh", "bcfarfield", "top", None),
            (ib, "klow", "bcsymmetryplane", "symm", None),
            (ib, "khigh", "bcsymmetryplane", "symm", None),
        ]
    return _write(blocks, bcs, fileName)


def sharp_fin(
    fileName,
    finAngle=15.0,
    xLE=-1.64,
    xEnd=0.25,
    yTop=0.15,
    zFar=0.26,
    ny=65,
    nz=89,
    dy0=2.0e-5,
    dz0=2.0e-5,
    Twall=300.0,
    xUp=0.02,
    ref=1.0,
):
    """Flat plate (y = 0) with a sharp fin whose apex is at x = 0.

    Freestream-aligned coordinates: x streamwise, y normal to the plate,
    z spanwise. The fin compression surface is z = x tan(finAngle) for
    x >= 0 and spans the full domain height. Upstream of the apex the
    z = 0 plane is a symmetry plane (undisturbed flow is along x).
    Blocks: 0 inviscid lead-in ahead of the plate, 1 plate up to the fin
    apex, 2 plate + fin.
    """
    ny = _n(ny - 1, ref)
    nz = _n(nz - 1, ref)
    xin = np.linspace(xLE - xUp, xLE, _n(9, ref))
    xall = fin_x_distribution(xLE, xEnd, ref)
    iApex = np.argmin(np.abs(xall))
    xall[iApex] = 0.0
    xs = [xin, xall[: iApex + 1], xall[iApex:]]
    y = tanh_one_sided(ny, yTop, dy0 / ref)
    tanA = np.tan(np.radians(finAngle))

    blocks = []
    for name, x in zip(["leadin", "plate", "fin"], xs):
        zb = np.maximum(x, 0.0) * tanA
        X = np.zeros((len(x), ny, nz, 3))
        for i, xi in enumerate(x):
            s = tanh_one_sided(nz, zFar - zb[i], dz0 / ref)
            X[i, :, :, 0] = xi
            X[i, :, :, 1] = y[:, None]
            X[i, :, :, 2] = zb[i] + s[None, :]
        blocks.append((name, X))

    bcs = [
        (0, "ilow", "bcinflowssupersonic", "inflow", None),
        (0, "jlow", "bcsymmetryplane", "symm_leadin", None),
        (1, "jlow", "bcwallviscousisothermal", "plate", _isothermal_dataset(Twall)),
        (2, "jlow", "bcwallviscousisothermal", "plate", _isothermal_dataset(Twall)),
        (0, "klow", "bcsymmetryplane", "symm_upstream", None),
        (1, "klow", "bcsymmetryplane", "symm_upstream", None),
        (2, "klow", "bcwallviscousisothermal", "fin", _isothermal_dataset(Twall)),
        (2, "ihigh", "bcoutflowsupersonic", "outflow", None),
    ]
    for ib in range(3):
        bcs += [(ib, "jhigh", "bcfarfield", "top", None), (ib, "khigh", "bcfarfield", "side", None)]
    return _write(blocks, bcs, fileName)


# ------------------------------------------------------------------------------
# Wideman et al. cylinder / offset flare (Mach 2.89)
# ------------------------------------------------------------------------------
R1 = 0.0254  # initial cylinder radius
R2 = 0.0635  # afterbody radius
ECC = 0.0127  # offset of the cone axis (towards phi = 0)
CONE = 20.0  # cone half-angle


def flare_radius(x, phi):
    """Body radius r_s(x, phi); x = 0 at the innermost (phi = 0) junction.
    phi = 0 is the side towards which the cone axis is displaced (+y)."""
    t = np.tan(np.radians(CONE))
    rho = (R1 - ECC) + x * t  # cone radius about its own axis
    arg = rho**2 - (ECC * np.sin(phi)) ** 2
    # Only the downstream nappe of the cone exists (rho > 0); upstream of its
    # apex the body is the initial cylinder.
    rc = np.where((arg > 0) & (rho > 0), ECC * np.cos(phi) + np.sqrt(np.maximum(arg, 0.0)), 0.0)
    return np.clip(np.maximum(R1, rc), R1, R2)


def junction_x(phi):
    """Axial location of the cylinder/cone junction, Eq. 30 of Edwards & Chandra."""
    return R1 * (np.sqrt(5.0 - 4.0 * np.cos(phi)) - 1.0) / (2.0 * np.tan(np.radians(CONE)))


def cylinder_flare(
    fileName,
    xLE=-0.90,
    xEnd=0.25,
    rOut=0.25,
    nr=81,
    nphi=49,
    dr0=3.0e-6,
    ref=1.0,
):
    """Cylinder with a sharp leading edge at xLE (the boundary layer develops
    along the cylinder), offset 20 deg flare and 12.7 cm afterbody.
    Half model: phi in [0, pi] with symmetry planes at phi = 0 and pi."""
    nr = _n(nr - 1, ref)
    nphi = _n(nphi - 1, ref)
    x = piecewise(
        [
            (xLE, -0.06, _n(41, ref), 2.0e-4 / ref, 1.0e-2 / ref),
            (-0.06, 0.13, _n(101, ref), 1.0e-2 / ref, 1.5e-3 / ref),
            (0.13, xEnd, _n(21, ref), 1.5e-3 / ref, 1.0e-2 / ref),
        ]
    )
    phi = np.linspace(0.0, np.pi, nphi)
    X = np.zeros((len(x), nr, nphi, 3))
    for k, p in enumerate(phi):
        rs = flare_radius(x, p)
        for i in range(len(x)):
            r = rs[i] + tanh_one_sided(nr, rOut - rs[i], dr0 / ref)
            X[i, :, k, 0] = x[i]
            X[i, :, k, 1] = r * np.cos(p)
            X[i, :, k, 2] = r * np.sin(p)
    blocks = [("body", X)]
    bcs = [
        (0, "ilow", "bcinflowssupersonic", "inflow", None),
        (0, "ihigh", "bcoutflowsupersonic", "outflow", None),
        (0, "jlow", "bcwallviscous", "body", None),
        (0, "jhigh", "bcfarfield", "far", None),
        (0, "klow", "bcsymmetryplane", "symm", None),
        (0, "khigh", "bcsymmetryplane", "symm", None),
    ]
    return _write(blocks, bcs, fileName)


if __name__ == "__main__":
    import sys, os

    out = sys.argv[1] if len(sys.argv) > 1 else "grids"
    os.makedirs(out, exist_ok=True)
    flat_plate_2d(os.path.join(out, "flat_plate_M82.cgns"))
