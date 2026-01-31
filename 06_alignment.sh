#!/usr/bin/env bash
#Requires: conda env "minimap2_env"
set -euo pipefail

MEDAKA_ASSEMBLY_PATH="$HOME/medaka_out/consensus.fasta"
REFERENCE="GCA_002507875.2_ASM250787v2_genomic.fna"
READS="SRR32410565.fastq.gz"

#Sanity Check
[[ -f "$READS" ]] || { echo "FASTQ not found: $READS" >&2; exit 1; }
[[ -f "$MEDAKA_ASSEMBLY_PATH" ]] || { echo "Medaka assembly not found: $MEDAKA_ASSEMBLY_PATH" >&2; exit 1; }
[[ -f "$REFERENCE" ]]      || { echo "Reference not found: $REFERENCE" >&2; exit 1; }


#Align raw reads to reference genome
minimap2 -a -x map-ont "$REFERENCE" "$READS" > aln.sam

#Convert to .BAM file, sort and index
samtools view -bS aln.sam -o aln.bam
samtools sort aln.bam -o alm.sorted.bam
samtools index alm.sorted.bam

#Align assembled genome to reference genome
minimap2 -a -x asm5 "$REFERENCE" "$MEDAKA_ASSEMBLY_PATH" > assembled_aln.sam

#Convert to .BAM file, sort and index
samtools view -bS assembled_aln.sam -o assembled_aln.bam
samtools sort assembled_aln.bam -o assembled_aln.sorted.bam
samtools index assembled_aln.sorted.bam
