#!/bin/bash
# Usage: ./run.sh NPROC script.py [args...]
# Runs a validation script with this ADflow build (not the pip-installed one).
cd "$(dirname "$0")"
export PYTHONPATH="$(cd .. && pwd)${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg            # avoid matplotlib probing a (possibly stale) display
export GFORTRAN_UNBUFFERED_ALL=y  # stream solver output to the log
NP=$1; shift
exec mpirun -np "$NP" python -u "$@"
