#!/bin/bash
#SBATCH --job-name=orca_worker
#SBATCH --partition=kemi1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=240:00:00
#SBATCH --output=logs/inp_%A_%a_orca.out
#SBATCH --error=logs/inp_%A_%a_orca.err
python ../../scripts/utils/submit_sp_orca.py --DB_PATH "$PWD/db/jobs.db" --INP_DIR "$PWD/INPUT_SP" --OUT_DIR "$PWD/Output_SP" --cpus 4 --mem 4000 --batch_size 200

