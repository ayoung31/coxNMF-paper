#!/bin/bash

#SBATCH --mem-per-cpu=1g
#SBATCH -n 50
#SBATCH -t 04:00:00
#SBATCH --output=out/mask_dijk.out

Rscript run_mask_Dijk.R