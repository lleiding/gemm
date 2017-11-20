#!/bin/sh

sbatch islandsim_hf.sh # high tol, full lnkg
sbatch islandsim_hn.sh # high tol, no lnkg
sbatch islandsim_lf.sh # low tol, full lnkg
sbatch islandsim_ln.sh # low tol, no lnkg
