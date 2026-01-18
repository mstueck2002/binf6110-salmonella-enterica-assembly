# BINF6110 Assignment 1: Using Bioinformatic Research Tools to Assemble the Bacterial Genome of *Salmonella enterica* 
## Introduction

Bacterial infections constitute a major public health burden and are responsible for significant levels of morbidity and mortality across human and animal populations. [1, 2] The diverse host invasion strategies employed by pathogenic bacteria, together with the rapid evolution of antibiotic resistance, highlight the need for a comprehensive understanding of the molecular and cellular mechanisms underlying bacterial virulence and host-pathogen interactions. Advances in high-throughput genome sequencing technologies, coupled with modern bioinformatics tools, have substantially facilitated bacterial genome analysis. [2] The availability of high-quality, accurately reconstructed genomes is essential for these analyses, emphasizing genome assembly as a fundamental step in bacterial genomics. [3]

*Salmonella enterica* is a Gram-negative, facultative intracellular bacterium with more than 2,600 distinct serovars. [4, 5] As a highly adaptable environmental microorganism, it occupies a wide range of reservoirs and represents a major pathogen of both humans and animals. Human infection can result in salmonellosis, which may be severe or, in some cases, fatal [6]. The *S. enterica* genome is approximately 4.8 Mb in size and is characterized by high gene density and the presence of diverse mobile genetic elements, including plasmids and prophages. [7, 8] Given its significant public health relevance and extensive serovar diversity, *S. enterica* is a well-suited model for bacterial genome assembly and comparative genomic analyses. Furthermore, numerous high-quality, well-annotated reference genomes are publicly available, allowing for effective validation of assembled genomes. [7]

Reconstructing a bacterial genome from sequencing reads requires technologies capable of achieving high sequence coverage, minimal fragmentation, and sufficient accuracy, while also maintaining short turnaround times to accommodate rapid bacterial evolution. [9, 10]  Oxford Nanopore Technologies (ONT) R10 sequencing has enabled the generation of near-complete microbial genomes with coverage exceeding 40-fold, in some cases without the need for additional short-read polishing. [11] Furthermore, ONT sequencing is rapid and efficient, supporting real-time base calling and analysis, which is particularly advantageous for time-sensitive applications such as outbreak response.[12] Despite these advantages, ONT-derived assemblies may exhibit systematic base-calling errors, particularly at methylated genomic positions, and long reads can span highly repetitive regions that complicate accurate assembly. [10, 13] To mitigate these limitations, post-assembly polishing strategies are commonly employed to improve consensus accuracy and produce high-quality microbial genome assemblies. [14] 

Before the advent of autonomous long-read assemblers, bacterial genome assembly was labour-intensive and limited in its ability to resolve repetitive genomic regions. [15] Flye was developed as an automated de novo assembler optimized for long-read sequencing data, enabling efficient and accurate genome reconstruction. The tool employs a repeat graph–based approach to resolve repetitive regions and generate contiguous genome assemblies [13]. Comparative studies have shown that Flye generates more accurate and contiguous assemblies than other state-of-the-art single-molecule sequencing assemblers, while also effectively reconstructing mosaic segmental duplication structures.¹³ Despite these advantages, base-calling errors introduced during Oxford Nanopore sequencing can propagate into the initial assembly. Post-assembly polishing tools address these residual errors, improving overall accuracy to levels comparable with short-read sequencing approaches. [14] Among Nanopore polishing methods, Medaka has been shown to produce high-quality, low-error genomes suitable for downstream analyses. [14]

Once a high-quality consensus genome has been assembled and polished, comparison to a well-annotated reference genome enables the identification and characterization of sequence variation.[16] Alignment to a reference strain facilitates the detection of single-nucleotide polymorphisms (SNPs) and insertions/deletions (indels), which are essential for understanding strain-level differences and evolutionary processes in bacterial genomes.[16] For reference-based mapping of long-read sequencing data, minimap2 is widely regarded as a benchmark method.[17] It employs a seed–chain–align strategy to map long sequencing reads to reference genomes, achieving high accuracy and computational efficiency relative to other domain-specific alignment tools.[18] Although minimap2 enables efficient alignment of large-scale genomic datasets, the resulting alignments can be challenging to interpret without visualization. The Integrative Genomics Viewer (IGV) addresses this by using efficient, multi-resolution file formats that support real-time exploration of large genomic datasets with minimal computational expense. [19] 

This study evaluates a computational pipeline for bacterial genome assembly using Oxford Nanopore R10 sequencing data to assemble and polish a Salmonella enterica genome, followed by reference-based comparison to identify and visualize sequence variation.

## Methods
### Sequencing Data Acquisition
Oxford Nanopore sequencing reads generated using R10 chemistry (expected accuracy Q20+, N50:5-15kb) were obtained in FASTQ format using the SRA Toolkit from NCBI (SRA accession: SRR32410565). Raw sequencing reads were used directly for downstream analyses without prior reference bias. Read quality and length distributions were assessed using standard Seqtk summary statistics to confirm suitability for long-read genome assembly. Ultra-short reads (<1000 bp), which provide limited structural information and can introduce noise into long-read assembly graphs, were filtered out with Seqtk.

### Genome Assembly
De novo genome assembly was performed using Flye (v2.9.2), an autonomous long-read assembler optimized for error-prone single-molecule sequencing data such as those generated by Oxford Nanopore sequencing. [13] Assembly was conducted using the --nano-hq parameter to specify high-quality Nanopore input reads, while all other parameters were retained at default settings to minimize manual bias. Flye employs a repeat graph-based assembly approach, enabling effective resolution of repetitive genomic regions and the generation of contiguous genome assemblies.

### Assembly Polishing
To improve consensus accuracy and mitigate base-calling errors associated with Oxford Nanopore sequencing, post-assembly polishing was performed using Medaka (v1.7.3) [14]. Medaka applies neural network–based models trained on Oxford Nanopore data to correct systematic sequencing errors. The appropriate Medaka model corresponding to R10 chemistry was selected (r1041_e82_400bps_sup_v5.0.0), and polishing was conducted using default settings to refine the assembled genome before downstream analyses.

### Reference Genome Selection
A high-quality, well-annotated Salmonella enterica reference genome (NCBI Assembly ASM250787v2) was retrieved from the NCBI Assembly database [7]. This reference genome was selected based on its completeness and annotation quality, providing a reliable framework for comparative analyses and variant detection.

### Read Alignment to Reference Genome
Polished sequencing reads were aligned to the reference genome using minimap2 (v2.26), a widely used and benchmarked aligner for long-read sequencing data [17]. Alignment was performed using the -ax map-ont preset, which is optimized for Oxford Nanopore reads and balances alignment accuracy with computational efficiency. Resulting alignments were stored in BAM format for downstream analysis.

### Variant Identification and Visualization
Aligned reads were used to identify sequence variation between the assembled genome and the reference genome. Variants of interest included single-nucleotide polymorphisms (SNPs) and insertions/deletions (indels), which provide insight into strain-level differences and evolutionary divergence.[16] Genomic alignments and detected variants were ​​visualized using IGV (v2.16.2) to enable interactive inspection of alignments and variants [19]

## Works Cited

1. Doron, S., & Gorbach, S. L. (2008). Bacterial Infections: Overview. International Encyclopedia of Public Health, 273–282. https://doi.org/10.1016/B978-012373960-5.00596-7

2. Donkor E. S. (2013). Sequencing of bacterial genomes: principles and insights into pathogenesis and development of antibiotics. Genes, 4(4), 556–572. https://doi.org/10.3390/genes4040556

3. Kumar, M. S., Krishna, M. B., Soman, K. P., Stanley, J., Pourmand, N., Suravajhala, P., & Babu, T. G. S. (2025). Benchmarking long-read assembly tools and preprocessing strategies for bacterial genomes: A case study on E. coli DH5α. Biotechnology reports (Amsterdam, Netherlands), 48, e00931. https://doi.org/10.1016/j.btre.2025.e00931 

4. Andino, A., Hanning, I., Salmonella enterica: Survival, Colonization, and Virulence Differences among Serovars, The Scientific World Journal, 2015, 520179, 16 pages, 2015. https://doi.org/10.1155/2015/520179 

5. Brown EW, Bell R, Zhang G, Timme R, Zheng J, Hammack TS, Allard MW.2021.Salmonella Genomics in Public Health and Food Safety. 9:eESP-0008-2020.https://doi.org/10.1128/ecosalplus.ESP-0008-2020 

6. Gourama, H. (2020). Foodborne Pathogens. In: Demirci, A., Feng, H., Krishnamurthy, K. (eds) Food Safety Engineering. Food Engineering Series. Springer, Cham. https://doi.org/10.1007/978-3-030-42660-6_2

7. National Center for Biotechnology Information (NCBI). Salmonella enterica genome assembly ASM250787v2. NCBI Assembly database. https://www.ncbi.nlm.nih.gov/assembly/ASM250787v2

8. Andrews, K., Landeryou, T., Sicheritz-Pontén, T., & Nale, J. Y. (2024). Diverse Prophage Elements of Salmonella enterica Serovars Show Potential Roles in Bacterial Pathogenicity. Cells, 13(6), 514. https://doi.org/10.3390/cells13060514

9. Wick, R. R., Judd, L. M., & Holt, K. E. (2023). Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. PLoS computational biology, 19(3), e1010905. https://doi.org/10.1371/journal.pcbi.1010905

10. Bogaerts B, Maex M, Commans F, Goeders N, Van den Bossche A, De Keersmaecker SCJ, Roosens NHC, Ceyssens P, Mattheus W, Vanneste K.2025.Oxford Nanopore Technologies R10 sequencing enables accurate cgMLST-based bacterial outbreak investigation of Neisseria meningitidis and Salmonella enterica when accounting for methylation-related errors. JClinMicrobiol63:e00410-25.https://doi.org/10.1128/jcm.00410-25

11. Sereika, M., Kirkegaard, R.H., Karst, S.M. et al. Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. Nat Methods 19, 823–826 (2022). https://doi.org/10.1038/s41592-022-01539-7

12. Kono N, Arakawa K. Nanopore sequencing: Review of potential applications in functional genomics. Develop Growth Differ. 2019; 61: 316–326. https://doi.org/10.1111/dgd.12608

13. Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8

14. Lee, J.Y., Kong, M., Oh, J. et al. Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis. Sci Rep 11, 20740 (2021). https://doi.org/10.1038/s41598-021-00178-w

15. Mihai Pop, Genome assembly reborn: recent computational challenges, Briefings in Bioinformatics, Volume 10, Issue 4, July 2009, Pages 354–366, https://doi.org/10.1093/bib/bbp026

16. Olson, N. D., Lund, S. P., Colman, R. E., Foster, J. T., Sahl, J. W., Schupp, J. M., Keim, P., Morrow, J. B., Salit, M. L., & Zook, J. M. (2015). Best practices for evaluating single nucleotide variant calling methods for Microbial Genomics. Frontiers in Genetics, 6. https://doi.org/10.3389/fgene.2015.00235 

17. Liyanage, K., Samarakoon, H., Parameswaran, S. et al. Efficient end-to-end long-read sequence mapping using minimap2-fpga integrated with hardware accelerated chaining. Sci Rep 13, 20174 (2023). https://doi.org/10.1038/s41598-023-47354-8

18. Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

19. Robinson, J. T., Thorvaldsdóttir, H., Winckler, W., Guttman, M., Lander, E. S., Getz, G., & Mesirov, J. P. (2011). Integrative genomics viewer. Nature biotechnology, 29(1), 24–26. https://doi.org/10.1038/nbt.1754



