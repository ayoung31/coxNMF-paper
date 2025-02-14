#!/bin/bash

#SBATCH --mem-per-cpu=1g
#SBATCH -n 30
#SBATCH -t 02:00:00
#SBATCH --output=out/mask_cptac.out

Rscript run_mask_CPTAC.R