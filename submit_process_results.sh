#!/bin/bash

DSTRING=20180402
BURNIN=8000

qsub -N proc qsub_process_results.sh --prefix=msp${DSTRING} --burnin=${BURNIN}
qsub -N proc_f qsub_process_results.sh --prefix=msp_f${DSTRING} --fix_allom2 --burnin=${BURNIN}
qsub -N proc_h qsub_process_results.sh --prefix=msp_h${DSTRING} --hetero --burnin=${BURNIN}
qsub -N proc_fh qsub_process_results.sh --prefix=msp_hf${DSTRING} --fix_allom2 --hetero --burnin=${BURNIN}
