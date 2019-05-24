#!/bin/bash
#SBATCH --job-name=circumbinary
#SBATCH --partition=q20
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=1-1024%40
#SBATCH --output=./slurm_logs/%A_%a.out

echo "========= Job started  at `date` =========="

echo "My jobid: $SLURM_JOB_ID"
echo "My array id: $SLURM_ARRAY_TASK_ID"

#mkdir slurm_logs 

models_per_job=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
s=$(echo "20000 + $models_per_job * ($SLURM_ARRAY_TASK_ID - 1)" | bc -l)

module load anaconda3

python3 sobol_dispatcher.py -s $s -N $models_per_job 

echo "========= Job finished at `date` =========="
#
