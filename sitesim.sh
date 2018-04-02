#!/bin/bash
#$ -N ss_1000
#$ -q "geo*"
#$ -pe omp 2
#$ -j y
#$ -t 1-54

/usr3/graduate/ashiklom/.singularity/Rscript scripts/multi_site_pda/04_site_simulation.R 1000 $SGE_TASK_ID
