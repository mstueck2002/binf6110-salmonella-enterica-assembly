#!/usr/bin/env bash
#Requires: conda env "qc_env"
set -euo pipefail

#Specify paths
ASSEMBLY="out_nanopore_flye/10-consensus/consensus.fasta"
REFERENCE="reference_fasta.fna"

#Sanity checks
[[ -f "$ASSEMBLY" ]]  || { echo "Assembly not found: $ASSEMBLY" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

#Quality check assembly – QUAST
quast.py "$ASSEMBLY" \
  -r "$REFERENCE" \
  -o out_quast

#Quality check assembly – BUSCO
busco -i "$ASSEMBLY" \
  -o out_busco \
  -m genome \
  --auto-lineage
