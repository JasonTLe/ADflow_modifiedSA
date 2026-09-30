"""
Tests for the SA-Edwards turbulence model and the SA compressibility
correction (useCompressibilitySA), see src/turbulence/saCorrections.F90.

These tests do not use reference files. Instead they check that
- the blockette (NK/ANK) and the regular residual routines agree,
- the forward mode AD agrees with finite differences,
- the forward, reverse and reverse-fast AD routines are consistent
  (dot-product tests),
- the new terms actually change the turbulence residual (and only it).
"""

# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_tutorial_wing

baseDir = os.path.dirname(os.path.abspath(__file__))

baseOptions = {
    "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
    "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
    "equationType": "RANS",
    "smoother": "DADI",
    "useNKsolver": True,
    "useBlockettes": True,
}

test_params = [
    {"name": "SA_Edwards", "options": {"turbulenceModel": "SA-Edwards"}},
    {"name": "SA_comp", "options": {"turbulenceModel": "SA", "useCompressibilitySA": True}},
    {"name": "SA_Edwards_comp", "options": {"turbulenceModel": "SA-Edwards", "useCompressibilitySA": True}},
]


def _dot(comm, a, b):
    return comm.allreduce(np.dot(np.ravel(a), np.ravel(b)))


def _norm(comm, a):
    return np.sqrt(_dot(comm, a, a))


@parameterized_class(test_params)
class TestSACorrections(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(baseOptions)
        options.update(self.options)

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)
        self.ap = copy.deepcopy(ap_tutorial_wing)
        self.CFDSolver.getResidual(self.ap)
        self.comm = self.CFDSolver.comm
        self.nw = 6

    def test_blockette_consistency(self):
        # The NK residual is evaluated with the blockette code; compare the
        # turbulence residual with the regular (non-blockette) block code,
        # which uses saCorr_block. (The mean flow residuals of the two paths
        # differ in ADflow independently of the turbulence model.)
        res1 = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)[:, 5]
        self.CFDSolver.setOption("useBlockettes", False)
        res2 = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)[:, 5]
        self.CFDSolver.setOption("useBlockettes", True)
        err = _norm(self.comm, res1 - res2) / _norm(self.comm, res1)
        self.assertLess(err, 1e-12)

    def test_correction_changes_turb_residual(self):
        # Scaling the compressibility correction must only change the
        # turbulence residual, and do so linearly in C5.
        if not self.CFDSolver.getOption("useCompressibilitySA"):
            self.skipTest("compressibility correction not used")
        res1 = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)
        self.CFDSolver.setOption("SAc5", 0.0)
        res0 = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)
        self.CFDSolver.setOption("SAc5", 7.0)
        res2 = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)
        self.CFDSolver.setOption("SAc5", 3.5)

        d1 = res1 - res0
        d2 = res2 - res0
        # Mean flow residuals are unaffected, up to the round-off level at
        # which repeated residual evaluations differ (the restart state is
        # converged, so the mean flow residual itself is tiny).
        self.CFDSolver.setOption("SAc5", 3.5)
        res1b = self.CFDSolver.getResidual(self.ap).reshape(-1, self.nw)
        noise = _norm(self.comm, res1b - res1)
        self.assertLessEqual(_norm(self.comm, d1[:, :5]), 10 * noise + 1e-14)
        # The turbulence residual changes...
        dNorm = _norm(self.comm, d1[:, 5])
        self.assertGreater(dNorm / _norm(self.comm, res0[:, 5]), 1e-8)
        # ...linearly in C5.
        self.assertLess(_norm(self.comm, d2[:, 5] - 2.0 * d1[:, 5]) / dNorm, 1e-8)
        # The destruction term is subtracted from the source and the
        # residual is -vol*source, so the change is non-negative.
        self.assertGreaterEqual(np.min(d1[:, 5]), -1e-10 * np.max(np.abs(d1[:, 5])))

    def _fdErr(self, resDot, h, comp, **kwargs):
        resDotFD = self.CFDSolver.computeJacobianVectorProductFwd(residualDeriv=True, mode="FD", h=h, **kwargs)
        a = resDot.reshape(-1, self.nw)[:, comp].flatten()
        b = resDotFD.reshape(-1, self.nw)[:, comp].flatten()
        return _norm(self.comm, a - b) / _norm(self.comm, a)

    def test_fwd_vs_fd(self):
        # Forward (one-sided) FD converges at first order towards the AD
        # result: check the error is small and drops ~10x per decade in h.
        wDot = self.CFDSolver.getStatePerturbation(321)
        resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        for comp, name in [(slice(0, 5), "flow"), (5, "turb")]:
            err1 = self._fdErr(resDot, 1e-10, comp, wDot=wDot)
            err2 = self._fdErr(resDot, 1e-11, comp, wDot=wDot)
            self.assertLess(err2, 5e-4, msg=name)
            self.assertGreater(err1 / err2, 5.0, msg=name)

        xVDot = self.CFDSolver.getSpatialPerturbation(314)
        resDot = self.CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)
        err1 = self._fdErr(resDot, 1e-8, 5, xVDot=xVDot)
        err2 = self._fdErr(resDot, 1e-9, 5, xVDot=xVDot)
        print("xVDot FD errors", err1, err2)
        self.assertLess(min(err1, err2), 1e-4)

    def test_dot_products(self):
        wDot = self.CFDSolver.getStatePerturbation(321)
        xVDot = self.CFDSolver.getSpatialPerturbation(314)
        resBar = self.CFDSolver.getStatePerturbation(123)

        resDot_w = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        resDot_x = self.CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)

        wBar, xVBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=resBar, wDeriv=True, xVDeriv=True)
        wBarFast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=resBar)

        fwd = _dot(self.comm, resDot_w, resBar)
        np.testing.assert_allclose(_dot(self.comm, wBar, wDot), fwd, rtol=1e-11)
        np.testing.assert_allclose(_dot(self.comm, wBarFast, wDot), fwd, rtol=1e-11)

        fwd = _dot(self.comm, resDot_x, resBar)
        np.testing.assert_allclose(_dot(self.comm, xVBar, xVDot), fwd, rtol=1e-11)


if __name__ == "__main__":
    unittest.main()
