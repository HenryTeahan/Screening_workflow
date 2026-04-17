#!/bin/bash
#SBATCH --job-name=input_worker
#SBATCH --partition=kemi1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=100:00:00
#SBATCH --array=1-10
#SBATCH --output=logs/inp_%A_%a.out
#SBATCH --error=logs/inp_%A_%a.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

micromamba run -n Complex python ../../scripts/utils/make_sp_input.py --DB_PATH "$PWD/db/jobs.db" --INP_DIR "$PWD/INPUT_SP" --EMBED_DIR "$PWD/Embed" --query "!r2scan-3c smd(thf) SP" --mem 4000 --cpus 4
