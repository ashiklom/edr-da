#!/bin/bash
#$ -l h_rt=07:00:00
#$ -j y
#$ -o logs_edr
#$ -t 1-50
#$ -v OMP_THREAD_LIMIT=1

/usr3/graduate/ashiklom/.singularity/Rscript scripts/ensemble/02_run_edr.R \
    --ens=$SGE_TASK_ID \
    $@
