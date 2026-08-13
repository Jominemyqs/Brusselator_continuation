#!/bin/bash
#SBATCH --job-name=target_mid
#SBATCH --time=144:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=4-5
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module load matlab

RUN_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$RUN_DIR"

export TASK="$SLURM_ARRAY_TASK_ID"

matlab-threaded -nodisplay -nosplash -r \
  "try, Main_target; catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"
