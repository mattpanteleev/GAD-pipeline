# GAD: Genome assembly, Annotation and Differential expression analysis
## What does this pipeline do?
This pipeline is specialised for assembling genome de novo for *S. cerevisiae* using short and long reads, annotating it using external evidence like RNA-seq data and analyzing the differences in gene expressions
## Pipeline Overview
This tool consists of 3 pipelines:
![](https://github.com/mattpanteleev/GAD-pipeline/blob/main/plots/whole%20pipeline.png)
- Genome Assembly
- Annotation
- Differential Expression Analysis
## Genome Assembly and Quality Control Pipeline
A genome assembly pipeline combining short-read (Illumina) and long-read (Nanopore) sequencing data with quality control steps and scaffolding.

The genome assembling pipeline:
![](https://github.com/mattpanteleev/GAD-pipeline/blob/main/plots/genome%20assembling.png)
### Trimming:
- Short-read with FastP
- Long-read with Chopper and Porechop
### Quality Control:
- Short-read QC with FastQC
- Long-read QC with NanoPlot
### Assembly:
- Initial assembly with Flye (long reads)
- Polishing with Hypo (short reads)
- Reference-based scaffolding with Ragtag

## Genome Annotation Pipeline
This Nextflow pipeline performs genome annotation using RNA-seq data, combining structural and functional annotation steps. It processes raw reads, aligns them to a reference genome, predicts gene structures, and integrates functional information from BLAST results and UniProt databases.

The genome annotation pipeline:
![](https://github.com/mattpanteleev/GAD-pipeline/blob/main/plots/genome%20annotation.png)

### Quality Control and Trimming:
- FastQC for initial read quality assessment
- Fastp for read trimming and filtering
- FastQC (post-trimming) for quality re-evaluation
### Alignment:
- STAR for genome indexing and RNA-seq read alignment
### Structural Annotation:
- GeMoMa for gene prediction using homology and RNA-seq evidence
### Functional Annotation:
- AGAT for merging functional annotations (BLAST results + UniProt)
- Gene naming and final GFF file generation

## RNA-Seq Analysis Pipeline
This Nextflow pipeline performs RNA-Seq data analysis, including quality control, read trimming, alignment, and gene-level quantification. After that R script performs
differential expression (DE) analysis and enrichment analysis.
**Pipeline Overview**

1. Quality control, alignment, gene quantification steps:
### Quality Control:
Initial and post-trimming QC reports using FastQC
### Read Trimming:
Adapter trimming, poly-G/X trimming, and quality filtering using fastp
### Alignment:
Genome indexing and alignment to a reference genome using STAR
### Quantification:
Gene-level read counting using featureCounts 

2. R-Based Analysis:
### Differential Expression: 
edgeR (exactTest for low-replicate data)
## Visualization: 
Volcano plots and GO enrichment (up/downregulated pathways)




# Prerequisites 
- Docker engine 1.10.x (or later)
- Nextflow 22.10 (or later)
- Conda 4.5 (or later)
# Quick Start
## Genome assembling
```
nextflow run genome_assembling.nf \
    --read1 reads_R1.fastq.gz \
    --read2 reads_R2.fastq.gz \
    --read_long nanopore_reads.fastq \
    --reference reference_genome.fasta \
    --outdir results
```
### Input Parameters
**Required:**
```
--read1              # Illumina R1 reads
--read2              # Illumina R2 reads 
--read_long          # Nanopore reads
--reference          # Reference genome
```
**Optional:**
```
--outdir             # Output directory (default: results)
--threads            # CPU threads (default: 4)
```

## Genome annotation

```
nextflow run genome_annotation.nf \
    --genome_fasta genome.fa \
    --read1 RNA_reads_R1.fastq.gz \
    --read2 RNA_reads_R2.fastq.gz \
    --ncbi_gff ref_genes.gff \
    --ncbi_genomic_fna ref_genome.fna \
    --uniprot_db uniprot_prot.fasta \
    --outdir results
```
### Input Parameters
**Required:**
```
--genome_fasta			# the novel genome from the previous step (FASTA)
--read1					# Forward reads (FASTQ)	
--read2					# Reverse reads (FASTQ)
--ncbi_gff				# Reference GFF file (for GeMoMa)
--ncbi_genomic_fna		# Reference genomic FASTA (for GeMoMa)
--uniprot_db			# UniProt database (for functional annotation)
```
**Optional:**
```
--outdir             # Output directory (default: results)
--threads            # CPU threads (default: 8)
```

## RNA-Seq Analysis
1. Run the Nextflow Pipeline
```
nextflow run get_count_matrix.nf \
  --genome_fasta genome.fa \
  --genome_gff annotation.gff \
  --read1 RNA_reads_R1.fastq.gz \
  --read2 RNA_reads_R2.fastq.gz \
  --outdir ./results
```
### Input Parameters
**Required:**
```
--genome_fasta			# the novel genome from the previous step (FASTA)
--genome_gff			# the novel genome annotation from the previous step (GFF)
--read1					# Forward reads (FASTQ)	
--read2					# Reverse reads (FASTQ)
```
**Optional:**
```
--outdir             # Output directory (default: results)
--threads            # CPU threads (default: 8)
```
2. Run the R Script:
Using output of featureCount for control and case conditions:
```
Rscript final_version_RNA_seq_edgeR.R
```
Ensure control_strain_gene_counts.txt and case_strain_gene_counts.txt are in the working directory or change names of variables in R file.


