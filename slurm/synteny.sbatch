#!/bin/bash
#SBATCH --job-name=synteny
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

module load StdEnv/2023
module load miniconda3

conda run -n syny bash scripts/07_synteny.sh
