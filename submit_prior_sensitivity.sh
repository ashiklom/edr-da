#!/bin/bash

DSTRING=20180402

qsub -N prior qsub_prior_sensivity.sh --prefix=msp${DSTRING}
qsub -N prior_f qsub_prior_sensivity.sh --prefix=msp_f${DSTRING} --fix_allom2
qsub -N prior_h qsub_prior_sensivity.sh --prefix=msp_h${DSTRING} --hetero
qsub -N prior_fh qsub_prior_sensivity.sh --prefix=msp_hf${DSTRING} --fix_allom2 --hetero
