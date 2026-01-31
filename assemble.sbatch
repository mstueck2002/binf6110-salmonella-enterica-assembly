#!/bin/bash
#SBATCH --job-name=flye_assembly
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

module load StdEnv/2023
source activate flye_env

bash scripts/02_assemble.sh
