#!/bin/bash -l
#PBS -l walltime=96:00:00  
#$ -S /bin/bash
#$ -cwd
#$ -N qsub_test
#$ -V

## Check modules are loaded
echo loading nescessary modules
module load openmpi/2.0.1-gnu540 hdf5/1.8.17-gcc540
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

R --vanilla < examples/inversion.R inversion.${log}.log 2> inversion.${log}.log
