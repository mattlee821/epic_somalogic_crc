#!/bin/bash

#SBATCH --job-name=analysis-combined
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /home/leem/001_projects/epic_somalogic_crc/

Rscript scripts/004_analysis/004_analysis.R
