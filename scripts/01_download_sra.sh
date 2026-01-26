

#!/usr/bin/env bash
set -euo pipefailh

#Use SRA tools to fetch sequence data
prefetch SRR32410565
fasterq-dump SRR32410565 --threads 12 

#Compress FASTQ file
gzip -f SRR32410565.fastq

#Verify file existence 
ls SRR32410565.fastq.gz

#View details 
zcat SRR32410565.fastq.gz 2>/dev/null | head
zcat SRR32410565.fastq.gz | wc -l
zcat SRR32410565.fastq.gz 2>/dev/null | awk 'NR%4==2 {print length($0); exit}'
seqtk fqchk SRR32410565.fastq.gz 2>/dev/null | head

