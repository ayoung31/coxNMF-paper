#!/bin/bash

#SBATCH --mem-per-cpu=200
#SBATCH -n 16
#SBATCH -t 00:50:00
#SBATCH --array=1-25
#SBATCH --output=out/info_std%a.out

sleep 10

R CMD BATCH "--args 16" ../run_stdNMF_sims.R out/stdnmf_sim${SLURM_ARRAY_TASK_ID}.Rout