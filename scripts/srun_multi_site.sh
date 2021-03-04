#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=s3673
#SBATCH --constraint='sky|hasw'
#SBATCH --output=logs/edr-da-revised-%j.log
#SBATCH --job-name=edr-da

source /home/ashiklom/.bash_functions
mod_r

module list

pwd; hostname; date

Rscript scripts/04-run-inversion.R $@
