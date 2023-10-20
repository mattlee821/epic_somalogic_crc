#!/bin/bash

#SBATCH --job-name=analysis-female-rectum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /home/leem/001_projects/epic_somalogic_crc/
  
Rscript scripts/003_analysis/rectum/001_analysis-female.R
