#!/bin/bash

## No resume
# sbatch --job-name=edr-plain --output=logs/edr-0.log srun_multi_site.sh
# sbatch --job-name=edr-h --output=logs/edr-h.log srun_multi_site.sh hetero
# sbatch --job-name=edr-ss --output=logs/edr-ss.log srun_multi_site.sh ss
# sbatch --job-name=edr-hss --output=logs/edr-hss.log srun_multi_site.sh hetero ss

# Resume
sbatch --job-name=edr-plain --output=logs/edr-0.log srun_multi_site.sh resume
sbatch --job-name=edr-h --output=logs/edr-h.log srun_multi_site.sh hetero resume
sbatch --job-name=edr-ss --output=logs/edr-ss.log srun_multi_site.sh ss resume
sbatch --job-name=edr-hss --output=logs/edr-hss.log srun_multi_site.sh hetero ss resume
