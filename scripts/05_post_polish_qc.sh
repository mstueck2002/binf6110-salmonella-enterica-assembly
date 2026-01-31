#!/usr/bin/env bash
#Requires: conda env "quast_env" and "busco_x86"
set -euo pipefail

MEDAKA_ASSEMBLY_PATH="$HOME/medaka_out/consensus.fasta"
REFERENCE="GCA_002507875.2_ASM250787v2_genomic.fna"

#Sanity check
[[ -f "$MEDAKA_ASSEMBLY_PATH" ]] || { echo "Medaka Assembly not found: $MEDAKA_ASSEMBLY_PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]]      || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

#Second quality check asssembly - QUAST
/opt/anaconda3/envs/quast_env/bin/quast.py "$MEDAKA_ASSEMBLY_PATH" \
  -r "$REFERENCE" \
  -o out_polished_quast

#Second quality check assembly - BUSCO
/opt/anaconda3/envs/busco_x86/bin/busco \
-i $MEDAKA_ASSEMBLY_PATH \
-o out_polished_busco \
-m genome \
--auto-lineage
