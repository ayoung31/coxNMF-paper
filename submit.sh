#!/bin/bash

#SBATCH --mem-per-cpu=1g
#SBATCH -n 80
#SBATCH -t 02:00:00
#SBATCH --array=1-25
#SBATCH --output=out/info%a.out

sleep 10

R CMD BATCH "--args 80" ../run_sims.R out/cv${SLURM_ARRAY_TASK_ID}.Rout