#!/bin/bash
#SBATCH --job-name=genome_polishing
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

module load StdEnv/2023
source activate medaka_env

bash scripts/04_polishing.sh
