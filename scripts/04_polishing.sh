#!/usr/bin/env bash
#Requires: conda env "medaka_env"
set -euo pipefail

ASSEMBLY_PATH="$HOME/out_nanopore_flye/assembly.fasta"
READS="SRR32410565.fastq.gz"
MODEL="r1041_e82_400bps_sup_v5.0.0"

#Sanity check
[[ -f "$READS" ]] || { echo "FASTQ not found: $READS" >&2; exit 1; }
[[ -f "$ASSEMBLY_PATH" ]] || { echo "Assembly not found: $ASSEMBLY_PATH" >&2; exit 1; }

medaka_consensus -m "$MODEL" \
  -i "$READS" \
  -d "$ASSEMBLY_PATH" \
  -o medaka_out \
  -t 12
