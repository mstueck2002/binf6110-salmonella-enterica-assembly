#!/usr/bin/env bash
# Requires: conda env "quast_env" and "busco_x86"
set -euo pipefail

MEDAKA_ASSEMBLY="medaka_out/consensus.fasta"
REFERENCE="reference_fasta.fna"

# Sanity checks
[[ -f "$MEDAKA_ASSEMBLY" ]] || { echo "Medaka assembly not found: $MEDAKA_ASSEMBLY" >&2; exit 1; }
[[ -f "$REFERENCE" ]]       || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

# Second quality check – QUAST
quast.py "$MEDAKA_ASSEMBLY" \
  -r "$REFERENCE" \
  -o out_polished_quast

# Second quality check – BUSCO
busco -i "$MEDAKA_ASSEMBLY" \
  -o out_polished_busco \
  -m genome \
  --auto-lineage
