#!/bin/bash
#SBATCH --job-name=spiral4
#SBATCH --time=120:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-3
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

module load matlab
cd ~/continuation_runs   # wherever Main_spiral.m lives

export TASK=$SLURM_ARRAY_TASK_ID

# Use threaded MATLAB if you want it (often faster for heavy linear algebra / FFTs)
matlab-threaded -nodisplay -nosplash -r "Main_target; exit;"
