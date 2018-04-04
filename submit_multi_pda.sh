#!/bin/bash

qsub -N msp_hf qsub_multi_pda.sh --hetero --fix_allom2 --prefix=msp_hf$(date +%Y%m%d)
qsub -N msp_f qsub_multi_pda.sh --fix_allom2 --prefix=msp_f$(date +%Y%m%d)
qsub -N msp_h qsub_multi_pda.sh --hetero --prefix=msp_h$(date +%Y%m%d)
qsub -N msp qsub_multi_pda.sh --prefix=msp$(date +%Y%m%d)
