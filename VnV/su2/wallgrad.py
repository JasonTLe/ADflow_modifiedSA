import sys, os, glob, numpy as np
sys.argv=[sys.argv[0]]
import importlib.util
from compare import read_vtk, U, rho, Tinf, mu, qref, run_dir
g, R = 1.4, 287.05
cp = g * R / (g - 1)
def muS(T): return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.55) / (T + 110.55)
Tw = 300.0
# ADflow
d = np.load(os.path.join(run_dir("adflow", "plate", "medium", "SA"), "cells.npz")); x, w = d["x"], d["w"]
xu = np.unique(np.round(x[:, 0], 8)); xa = xu[np.argmin(np.abs(xu + 0.11))]
m = (np.abs(x[:, 0] - xa) < 1e-7) & (x[:, 2] < x[:, 2].min() + 1e-9)
o = np.argsort(x[m, 1]); yA = x[m, 1][o]; wA = w[m][o]
pA = 0.4 * (wA[:, 4] - 0.5 * wA[:, 0] * (wA[:, 1] ** 2 + wA[:, 2] ** 2)); TA = pA / wA[:, 0] * Tinf  # ADflow: p_inf/rho_inf = 1 (nondimensional)
kw = muS(Tw) * cp / 0.72
print("ADflow: y1=%.3e T1=%.2f K  -> k dT/dy (1-sided, T at wall=300) = %.4f W/cm2" % (yA[0], TA[0], kw * (TA[0] - Tw) / yA[0] / 1e4))
# 2nd-order one-sided from first two cell centres assuming T(0)=Tw: fit T = Tw + a y + b y^2
A = np.array([[yA[0], yA[0]**2], [yA[1], yA[1]**2]]); a, b = np.linalg.solve(A, [TA[0]-Tw, TA[1]-Tw])
print("ADflow quadratic fit: k dT/dy = %.4f W/cm2  (T1=%.2f, T2=%.2f at y=%.2e, %.2e)" % (kw * a / 1e4, TA[0], TA[1], yA[0], yA[1]))
# SU2
snap = sorted(glob.glob(os.path.join(run_dir("su2", "plate", "medium", "SA"), "flow*_*.vtk")), key=os.path.getmtime)[-1]; P, s = read_vtk(snap)
xs = np.unique(np.round(P[:, 0], 9)); xb = xs[np.argmin(np.abs(xs - xa))]
mm = np.abs(P[:, 0] - xb) < 1e-9; oo = np.argsort(P[mm, 1]); yS = P[mm, 1][oo]; TS = s["Temperature"][mm][oo]
print("SU2: node0 y=%.1e T=%.2f; node1 y=%.3e T=%.2f; node2 y=%.3e T=%.2f" % (yS[0], TS[0], yS[1], TS[1], yS[2], TS[2]))
A = np.array([[yS[1], yS[1]**2], [yS[2], yS[2]**2]]); a2, b2 = np.linalg.solve(A, [TS[1]-TS[0], TS[2]-TS[0]])
print("SU2 quadratic fit: k dT/dy = %.4f W/cm2;  SU2 reported Heat_Flux at wall node = %.4f W/cm2" % (kw * a2 / 1e4, abs(s["Heat_Flux"][mm][oo][0]) / 1e4))
print("SU2 laminar viscosity at wall node: %.4e (Sutherland(300K) = %.4e)" % (s["Laminar_Viscosity"][mm][oo][0], muS(300.0)))
