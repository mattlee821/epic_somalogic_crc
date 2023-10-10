#!/bin/bash

#SBATCH --job-name=metaboprep-epic-somalogic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20-1:0:00
#SBATCH --mem=100000M

# SCRIPT=/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/scripts/002_metaboprep/run_metaboprep_pipeline.R
# DATA=/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/scripts/002_metaboprep/
# IN=/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/analysis/001_phenofile/
# OUT=/Users/leem/Library/CloudStorage/OneDrive-IARC/001_projects/epic_metabolomics_proteomics/analysis/002_metaboprep/
SCRIPT=/home/leem//001_projects/epic_somalogic_crc/scripts/002_metaboprep/run_metaboprep_pipeline.R
DATA=/home/leem/001_projects/epic_somalogic_crc/scripts/002_metaboprep/
IN=/home/leem/001_projects/epic_somalogic_crc/analysis/001_phenofile/
OUT=/home/leem/001_projects/epic_somalogic_crc/analysis/002_metaboprep/

# if running on HPC likely to get an error at the end when generating a report. to make report, copy newly made folder locally and run:
# output_dir_path = paste0("FULL/PATH/TO/DIR/OF/CHOICE/")
# rdfile = paste0(output_dir_path, "ReportData.Rdata")
# generate_report( full_path_2_Rdatafile = rdfile, dir_4_report = output_dir_path )
# when finished, copy over the new html report from metaboprep directory AND in the script location the new .md file and the figures/

# cancer
cd ${DATA}/cancer
Rscript ${SCRIPT} ${DATA}cancer/parameter_file_female.txt
mv ${IN}cancer/metaboprep* ${OUT}cancer/female/
mv ${DATA}cancer/figure/ ${OUT}cancer/female/
mv ${DATA}cancer/metaboprep* ${OUT}cancer/female/
