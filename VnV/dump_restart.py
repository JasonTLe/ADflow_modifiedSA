"""Write cells.npz (cell centres + states) of an existing solution.
Usage: ./run.sh N dump_restart.py {fin,flare} GRID.cgns VOLFILE OUT.npz"""
import sys, os
from mpi4py import MPI
from adflow import ADFLOW
from common import HERE, kussoy_ap, wideman_ap, solver_options, dump_cells
case, grid, vol, outf = sys.argv[1:5]
ap = kussoy_ap(case) if case in ("fin", "plate") else wideman_ap("flare")
o = solver_options(os.path.join(HERE, "grids", grid), os.path.dirname(os.path.abspath(outf)), "SA",
                   restartFile=vol, writeVolumeSolution=False, writeSurfaceSolution=False, writeTecplotSurfaceSolution=False)
S = ADFLOW(options=o, comm=MPI.COMM_WORLD)
S.getResidual(ap)
dump_cells(S, outf)
