#!/bin/bash

#SBATCH --mem=500
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH --array=1-100
#SBATCH --output=out/sim%a_nu10.out

R CMD BATCH "--args 10" sim_dat.R out/sim${SLURM_ARRAY_TASK_ID}_nu10.Rout