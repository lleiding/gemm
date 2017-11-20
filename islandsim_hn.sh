#!/bin/bash
#SBATCH -n 12 # Number of cores
#SBATCH --mem-per-cpu=8G # Memory pool per core/task (see also --mem for all cores)
#SBATCH -o %N_%A_%a.out # File to which STDOUT will be written %a: comination of array and job %j
#SBATCH -e %N_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ludwig.leidinger@uni-wuerzburg.de # Email to which notifications will be sent

/home/s216849/builds/julia-903644385b/bin/julia -p 4 main.jl -m mapfile1,mapfile"${SLURM_ARRAY_TASK_ID}" -s 1 -t high -l none -e low -c low
