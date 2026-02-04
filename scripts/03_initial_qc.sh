#!/usr/bin/env bash
#Requires: conda env "quast_env" and "busco_x86"
set -euo pipefail

#Specify paths
ASSEMBLY_PATH="$HOME/out_nanopore_flye/assembly.fasta"
REFERENCE="reference_fasta.fna"

#Sanity check
[[ -f "$ASSEMBLY_PATH" ]] || { echo "Assembly not found: $ASSEMBLY_PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]]      || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

#Quality check asssembly - QUAST
/opt/anaconda3/envs/quast_env/bin/quast.py "$ASSEMBLY_PATH" \
  -r "$REFERENCE" \
  -o out_quast

#Quality check assembly - BUSCO
/opt/anaconda3/envs/busco_x86/bin/busco \
-i $ASSEMBLY_PATH \
-o out_busco \
-m genome \
--auto-lineage

