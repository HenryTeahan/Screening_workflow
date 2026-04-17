#!/bin/bash
#SBATCH --job-name=embed_worker
#SBATCH --partition=kemi1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4000
#SBATCH --time=100:00:00
#SBATCH --array=1-50
#SBATCH --output=logs/embed_%A_%a.out
#SBATCH --error=logs/embed_%A_%a.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs
python ../../scripts/utils/embed_worker.py --DB_PATH "$PWD/db/jobs.db" --EMBED_DIR "$PWD/Embed" --N_tries 25

