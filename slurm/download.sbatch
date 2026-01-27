#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out

module load StdEnv/2023
source activate sra_env

bash scripts/01_download_sra.sh
