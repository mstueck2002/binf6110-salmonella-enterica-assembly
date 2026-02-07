#!/bin/bash
#SBATCH --job-name=post_polish_quality_check
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

module load StdEnv/2023
module load miniconda3

#Run second QC (QUAST + BUSCO)
conda run -n qc_env bash scripts/05_post_polish_qc.sh
