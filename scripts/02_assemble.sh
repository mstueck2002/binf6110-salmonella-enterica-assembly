#!/usr/bin/env bash
#Requires: conda env "flye_env"
set -euo pipefail

#Run Flye assembler
flye --nano-hq SRR32410565.fastq.gz --out-dir out_nanopore_flye --threads 12

#Sanity check
[[ -f SRR32410565.fastq.gz ]] || { echo "FASTQ not found"; exit 1; }
