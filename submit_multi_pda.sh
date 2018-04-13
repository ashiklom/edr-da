#!/bin/bash

DSTRING=20180402

qsub -N msp_hf qsub_multi_pda.sh --hetero --fix_allom2 --prefix=msp_hf$DSTRING
qsub -N msp_f qsub_multi_pda.sh --fix_allom2 --prefix=msp_f$DSTRING
qsub -N msp_h qsub_multi_pda.sh --hetero --prefix=msp_h$DSTRING
qsub -N msp qsub_multi_pda.sh --prefix=msp$DSTRING
