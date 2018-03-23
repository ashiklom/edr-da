#!/bin/bash
#$ -pe omp 16
#$ -j y
#$ -o logs
#$ -l h_rt=96:00:00

mkdir -p logs
#/usr3/graduate/ashiklom/.singularity/Rscript scripts/multi_site_pda/01_run_ed2.R
/usr3/graduate/ashiklom/.singularity/Rscript scripts/multi_site_pda/02_run_edr_pda.R 16
