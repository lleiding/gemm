#!/bin/sh

sbatch --array=3-6 islandsim_hf.sh
sbatch --array=3-6 islandsim_hn.sh
sbatch --array=3-6 islandsim_lf.sh
sbatch --array=3-6 islandsim_ln.sh
