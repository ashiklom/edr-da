#!/bin/bash
#$ -j y
#$ -o logs
#$ -pe omp 16
#$ -N dlsite

/usr3/graduate/ashiklom/.singularity/Rscript scripts/download_site_narr.R
