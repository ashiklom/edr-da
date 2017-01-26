#!/bin/bash -l
#PBS -l walltime=96:00:00  
#$ -S /bin/bash
#$ -cwd
#$ -N qsub_test
#$ -V
#$ -M sserbin@bnl.gov

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

## Run R PDA script
#cd ../
#R --vanilla < inversion_ps.prior.R > inversion_ps.prior.log 2> inversion_ps.prior.log
R --vanilla < inversion_ps.prior.R > inversion_ps.prior.1.log 2> inversion_ps.prior.1.log
#R --vanilla < inversion_ps.prior_fb.R > inversion_ps.prior_fb.log 2> inversion_ps.prior_fb.log
##


