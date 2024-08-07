#!/bin/bash

#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 00:04:00
#SBATCH --array=1-100
#SBATCH --output=out/res%a.out

sleep 10

R CMD BATCH ../full_models.R out/res${SLURM_ARRAY_TASK_ID}.Rout