#!/bin/bash -l
#PBS -l walltime=96:00:00  
#$ -S /bin/bash
#$ -cwd
#$ -N qsub_test
#$ -V

## Check modules are loaded
echo Checking loaded modules
module list
##

## Identify working environment
echo Starting directory: $PWD
cd $PBS_O_WORKDIR
echo Working directory: $PWD
echo $PBS_SERVER
echo $PBS_O_HOST 
##

## Log number
echo Log number: $log
##

## Run R PDA script

# original inversion_ps.prior.R
#R --vanilla < inversion_ps.prior.R > inversion_ps.prior.${log}.log 2> inversion_ps.prior.${log}.log

# updated multi-variate prior version
R --vanilla < inversion_ps.singlePFT.multivariate.prior.R inversion_ps.singlePFT.multivariate.prior.${log}.log 2> inversion_ps.singlePFT.multivariate.prior.${log}.log
##


