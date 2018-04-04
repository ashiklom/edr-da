#!/bin/bash
#$ -pe omp 2
#$ -j y
#$ -o logs_prior_sens/
#$ -l h_rt=120:00:00
#$ -t 1-54

mkdir -p logs_prior_sens

/usr3/graduate/ashiklom/.singularity/Rscript \
    scripts/prior_sensitivity/01.run_prior_simulation.R \
    --geo \
    --nsim=1000 \
    --nsite=${SGE_TASK_ID} \
    $@
