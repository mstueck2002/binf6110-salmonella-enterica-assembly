#!/usr/bin/env bash
# Requires: conda env "medaka_env"
set -euo pipefail

ASSEMBLY="out_nanopore_flye/10-consensus/consensus.fasta"
FASTQ="SRR32410565.fastq.gz"
MODEL="r1041_e82_400bps_sup_v5.0.0"
THREADS=10

# Sanity checks
[[ -f "$FASTQ" ]]    || { echo "FASTQ not found: $FASTQ" >&2; exit 1; }
[[ -f "$ASSEMBLY" ]] || { echo "Assembly not found: $ASSEMBLY" >&2; exit 1; }

# Run Medaka polishing
medaka_consensus -m "$MODEL" \
  -i "$FASTQ" \
  -d "$ASSEMBLY" \
  -o medaka_out \
  -t "$THREADS"
