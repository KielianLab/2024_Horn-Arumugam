#!/bin/bash
#SBATCH --time=167:00:00
#SBATCH --mem-per-cpu=500000
#SBATCH --job-name=kielian_OrthoInt_DE

module load R/4.0
R CMD BATCH HCC_2022_ortho_integration_DE.R
