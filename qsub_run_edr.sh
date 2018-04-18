#!/bin/bash
#$ -q "geo*"
#$ -pe omp 16
#$ -l h_rt=04:00:00
#$ -j y
#$ -o logs_edr

/usr3/graduate/ashiklom/.singularity/Rscript scripts/ensemble/02_run_edr.R \
  --ncores=16 \
  $@
