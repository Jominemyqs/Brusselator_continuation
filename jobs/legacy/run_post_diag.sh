#!/bin/bash
#SBATCH --job-name=postdiag
#SBATCH --output=logs/postdiag_%j.out
#SBATCH --error=logs/postdiag_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

mkdir -p logs diag_outputs

module load matlab
matlab -nodisplay -r "addpath(pwd); post_diagnostics_fast_target('target_right_up.mat', 'diag_outputs', 20); exit"
