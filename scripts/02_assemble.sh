#!/usr/bin/env bash
#Requires: conda env "flye_env"
set -euo pipefail

READS="SRR32410565.fastq.gz"

#Run Flye assembler
flye --nano-hq "$READS" --out-dir out_nanopore_flye --threads 12

#Sanity check
[[ -f "$READS" ]] || { echo "FASTQ not found"; exit 1; }
