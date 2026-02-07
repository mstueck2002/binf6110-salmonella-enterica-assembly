#!/usr/bin/env bash
#Requires: conda env "syny_env"
set -euo pipefail

MEDAKA_ASSEMBLY="medaka_out/consensus.fasta"
REFERENCE="reference_fasta.fna"

#Sanity checks
[[ -f "$MEDAKA_ASSEMBLY" ]] || { echo "Medaka assembly not found: $MEDAKA_ASSEMBLY" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

#Annotate the Medaka assembly
prokka \
  --outdir prokka_medaka \
  --prefix medaka \
  --kingdom Bacteria \
  --genus Salmonella \
  --species enterica \
  --usegenus \
  "$MEDAKA_ASSEMBLY"

#Annotate the reference genome - file structure similarity
prokka \
  --outdir prokka_reference \
  --prefix reference \
  --kingdom Bacteria \
  --genus Salmonella \
  --species enterica \
  --usegenus \
  "$REFERENCE"

#Synteny analysis between reference and Medaka assembly
run_syny.pl \
  -a prokka_reference/reference.gbk \
  -a prokka_medaka/medaka.gbk \
  -o syny_out
