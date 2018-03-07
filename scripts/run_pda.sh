#!/bin/bash
#$ -pe omp 4
#$ -j y
#$ -o logs/
#$ -N edr_pda
#$ -t 1-13
#$ -l h_rt=96:00:00

echo "==============================================="
echo "Current working directory: $PWD"
echo "===================== Starting PDA ============"
/usr3/graduate/ashiklom/.singularity/Rscript scripts/pda.R $SGE_TASK_ID
echo "===================== Done! ==================="
