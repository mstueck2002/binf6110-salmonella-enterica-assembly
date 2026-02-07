#!/usr/bin/env bash
#Requires: conda env "minimap2_env"
set -euo pipefail

MEDAKA_ASSEMBLY="medaka_out/consensus.fasta"
REFERENCE="reference_fasta.fna"
FASTQ="SRR32410565.fastq.gz"

#Sanity checks
[[ -f "$FASTQ" ]] || { echo "FASTQ not found: $FASTQ" >&2; exit 1; }
[[ -f "$MEDAKA_ASSEMBLY" ]] || { echo "Medaka assembly not found: $MEDAKA_ASSEMBLY" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "Reference not found: $REFERENCE" >&2; exit 1; }

#Align raw reads to reference genome
minimap2 -a -x map-ont "$REFERENCE" "$FASTQ" > aln.sam

#Convert to BAM, sort, and index
samtools view -bS aln.sam -o aln.bam
samtools sort aln.bam -o aln.sorted.bam
samtools index aln.sorted.bam

#Variant calling
bcftools mpileup \
  -f "$REFERENCE" \
  -Ou \
  -Q 7 \
  -q 0 \
  -d 10000 \
  -a FORMAT/DP,FORMAT/AD \
  aln.sorted.bam | \
bcftools call \
  -m \
  --ploidy 1 \
  -Ov \
  -o raw_calls.vcf

#Remove low-confidence calls
bcftools filter \
  -e 'QUAL<20 || FORMAT/DP<10' \
  raw_calls.vcf \
  -Ov -o filtered.vcf

#Keep only variant sites
bcftools view \
  -v snps,indels \
  filtered.vcf \
  -Ov -o variants_only.vcf

#Summary statistics
bcftools stats variants_only.vcf | grep -E 'SN|INDEL'
bcftools stats variants_only.vcf | grep 'TSTV'

#Align assembled genome to reference
minimap2 -a -x asm5 "$REFERENCE" "$MEDAKA_ASSEMBLY" > assembled_aln.sam

#Convert to BAM, sort, and index
samtools view -bS assembled_aln.sam -o assembled_aln.bam
samtools sort assembled_aln.bam -o assembled_aln.sorted.bam
samtools index assembled_aln.sorted.bam

