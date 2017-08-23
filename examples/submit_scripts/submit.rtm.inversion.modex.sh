#!/bin/bash -l
#PBS -l walltime=96:00:00  
#$ -S /bin/bash
#$ -cwd
#$ -N qsub_test
#$ -V

## Check modules are loaded
echo loading nescessary modules
module load openmpi/2.1.1-gnu540 hdf5/1.8.19-gcc540
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

#R --vanilla < examples/inversion.R inversion.${log}.log 2> inversion.${log}.log

#R --vanilla < examples/inversion_single.R inversion.${log}.log 2> inversion.${log}.log


R --vanilla < examples/inversion_bayestools.R inversion_bayestools.${log}.log 2> inversion_bayestools.${log}.log
