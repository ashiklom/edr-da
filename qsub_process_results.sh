#!/bin/bash
#$ -j y
#$ -o logs_process_results/

/usr3/graduate/ashiklom/.singularity/Rscript \
    scripts/multi_site_pda/03_process_results.R \
    $@
