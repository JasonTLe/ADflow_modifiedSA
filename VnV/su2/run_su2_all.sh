#!/bin/bash
# SU2 flat-plate runs: stage 1 (first order) then stage 2 (second order).
cd "$(dirname "$0")/../results/su2/plate"
SU2=/home/jason/packages/SU2/install/bin/SU2_CFD
for t in coarse medium fine; do for m in SA SA-Edwards; do
  d=$t/$m
  [ -f $d/stage2.log ] && grep -q "^EXIT 0" $d/stage2.log && continue
  (cd $d && timeout 30000 mpirun -np ${NP:-4} $SU2 stage1.cfg > stage1.log 2>&1; echo EXIT $? >> stage1.log
   mkdir -p stage1 && cp restart_flow.csv flow.vtk stage1/ 2>/dev/null
   timeout 30000 mpirun -np ${NP:-4} $SU2 stage2.cfg > stage2.log 2>&1; echo EXIT $? >> stage2.log)
done; done
