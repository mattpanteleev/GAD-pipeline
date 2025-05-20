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

# Prerequisites 
- Docker engine 1.10.x (or later)
- Nextflow 22.10 (or later)
- Conda 4.5 (or later, optional)
# Quick Start
```
nextflow run genome_assembling.nf \
    --read1 reads_R1.fastq.gz \
    --read2 reads_R2.fastq.gz \
    --read_long nanopore_reads.fastq \
    --reference reference_genome.fasta \
    --outdir results
```
## Input Parameters
### Required:
```
--read1              # Illumina R1 reads
--read2              # Illumina R2 reads 
--read_long          # Nanopore reads
--reference          # Reference genome
```
### Optional:
```
--outdir             # Output directory (default: results)
--threads            # CPU threads (default: 4)
```
