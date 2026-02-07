#!/usr/bin/env bash
#Requires: conda env "flye_env"
set -euo pipefail

FASTQ="SRR32410565.fastq.gz"

#Sanity check
[[ -f "$FASTQ" ]] || { echo "FASTQ not found"; exit 1; }

#Run Flye assembler
flye --nano-hq "$FASTQ" --out-dir out_nanopore_flye --threads 12

