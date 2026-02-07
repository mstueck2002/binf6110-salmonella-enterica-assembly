#!/bin/bash
#SBATCH --job-name=initial_quality_check
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

module load StdEnv/2023
module load miniconda3

conda run -n qc_env bash scripts/03_initial_qc.sh
