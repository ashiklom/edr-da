#!/bin/bash
#$ -q "geo*"
#$ -j y
#$ -o logs_ed
#$ -t 1-200
#$ -N BH07_site_1-25669

printf -v ENSDIR "ens_%03d" $SGE_TASK_ID

/usr3/graduate/ashiklom/.singularity/sexec /projectnb/dietzelab/ashiklom/ED2/ED/build/ed_2.1-dbg -f /projectnb/dietzelab/ashiklom/edr-da/ensemble_outputs/msp_hf20180402/BH07_site_1-25669/$ENSDIR/ED2IN
