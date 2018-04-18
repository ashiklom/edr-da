#!/bin/bash
#$ -l h_rt=04:00:00
#$ -j y
#$ -o logs_edr
#$ -t 1-50

/usr3/graduate/ashiklom/.singularity/Rscript scripts/ensemble/02_run_edr.R \
    --ens=$SGE_TASK_ID \
    $@
