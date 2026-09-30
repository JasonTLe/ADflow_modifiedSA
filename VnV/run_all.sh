#!/bin/bash
# Regenerate all V&V results (grids, ADflow runs, SU2 runs, post-processing).
# This is the recipe used for the reported results; it takes many hours on a
# 16-core machine. Run sections individually as needed:
#   ./run_all.sh grids | plate | fin | flare | su2 | post
set -e
cd "$(dirname "$0")"
export MPLBACKEND=Agg
CPL='{"ANKCoupledSwitchTol": 1e-3, "L2Convergence": 1e-10}'   # fin: early coupled ANK
DEC='{"ANKCoupledSwitchTol": 1e-16, "L2Convergence": 1e-10}'  # coarse flare: decoupled ANK

grids() {
  python -c "
import numpy as np, meshes
r = np.sqrt(2)
meshes.sharp_fin('grids/sharp_fin_M82.cgns'); meshes.sharp_fin('grids/sharp_fin_M82_refc.cgns', ref=1/r)
meshes.cylinder_flare('grids/cylinder_flare_M289.cgns'); meshes.cylinder_flare('grids/cylinder_flare_M289_refc.cgns', ref=1/r)
for tag, ref in [('', 1.0), ('_refc', 0.5), ('_reff', 2.0)]:
    meshes.flat_plate_2d(f'grids/flat_plate_M82_LE164{tag}.cgns', xLE=-1.64, ref=ref)"
}

plate() {  # 2-D Mach 8.2 plate: BL initial condition + 1500 RK iterations + ANK (defaults)
  for m in SA SA-Edwards SA-comp; do ./run.sh 4 run_flatplate.py $m flat_plate_M82_LE164.cgns > logs/fp164_$m.log 2>&1; done
  for t in refc reff; do for m in SA SA-Edwards; do  # grid study only (coarse / fine)
    ./run.sh 4 run_flatplate.py $m flat_plate_M82_LE164_$t.cgns > logs/fp164_${m}_$t.log 2>&1
  done; done
  # Sensitivity checks cited in the README (-> results/adflow/plate/sensitivity/SA<RUN_TAG>)
  RUN_TAG=_vort ADFLOW_EXTRA='{"turbulenceProduction": "vorticity"}' ./run.sh 4 run_flatplate.py SA flat_plate_M82_LE164.cgns > logs/fp164_SA_vort.log 2>&1
  RUN_TAG=_nu3 ADFLOW_EXTRA='{"eddyVisInfRatio": 0.21044}' ./run.sh 4 run_flatplate.py SA flat_plate_M82_LE164.cgns > logs/fp164_SA_nu3.log 2>&1
  python dump_restart.py plate flat_plate_M82_LE164.cgns results/adflow/plate/medium/SA/flatplate_000_vol.cgns results/adflow/plate/medium/SA/cells.npz
}

fin() {  # Mach 8.2 sharp fin: BL initial condition + 1500 RK + ANK
  ADFLOW_EXTRA="$CPL" ./run.sh 16 run_case.py fin SA > logs/fin_SA.log 2>&1
  ADFLOW_EXTRA="$CPL" ./run.sh 16 run_case.py fin SA-Edwards results/adflow/fin/medium/SA/fin_000_vol.cgns > logs/fin_SA-Edwards.log 2>&1
  for m in SA SA-Edwards; do  # coarse grid, cold start (restarts across models can stall)
    GRID=sharp_fin_M82_refc.cgns ADFLOW_EXTRA="$CPL" ./run.sh 8 run_case.py fin $m > logs/fin_${m}_refc.log 2>&1
  done
  for l in medium coarse; do for m in SA SA-Edwards; do python postprocess.py fin $l $m > /dev/null; done; done
}

flare() {  # Mach 2.89 cylinder / flare: BL initial condition + ANK
  for m in SA SA-Edwards SA-comp; do ./run.sh 12 run_case.py flare $m > logs/flare_$m.log 2>&1; done
  for m in SA SA-Edwards; do  # coarse grid: the coupled ANK solver stalls; keep it decoupled
    GRID=cylinder_flare_M289_refc.cgns STAGE1_RK=1500 ADFLOW_EXTRA="$DEC" ./run.sh 6 run_case.py flare $m > logs/flare_${m}_refc.log 2>&1
  done
  for m in SA SA-Edwards SA-comp; do python postprocess.py flare medium $m > /dev/null; done
  for m in SA SA-Edwards; do python postprocess.py flare coarse $m > /dev/null; done
}

su2() {  # SU2 v8.5.0 flat plate on the same grids (SU2 built in ~/packages/SU2)
  (cd su2 && for t in coarse medium fine; do for m in SA SA-Edwards; do python make_su2_case.py "$t" $m; done; done && ./run_su2_all.sh)
}

post() {
  python summarize.py > /dev/null
  python grid_study.py > /dev/null
  python plot_results.py
  (cd su2 && python compare.py > /dev/null && python profiles.py > /dev/null)
}

mkdir -p grids logs results/figures
steps=("$@")
[ ${#steps[@]} -eq 0 ] && steps=(grids plate fin flare su2 post)
for s in "${steps[@]}"; do $s; done
