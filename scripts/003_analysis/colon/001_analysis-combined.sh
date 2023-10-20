#!/bin/bash

#SBATCH --job-name=analysis-combined-colon
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

cd /home/leem/001_projects/epic_somalogic_crc/
  
Rscript scripts/003_analysis/colon/001_analysis-combined.R