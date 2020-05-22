#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --account=s3673

source /home/ashiklom/.bash_functions
mod_r

module list

pwd; hostname; date

Rscript scripts/inversion_r_multi_site.R $@
