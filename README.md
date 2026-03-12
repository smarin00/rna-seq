# RNA-seq: Ribosome Profiling

**Course:** Bioinformatics Master RNA-seq — MSc Bioinformatics & Computational Biology, Universität Bern
**Author:** Santiago Marin Martinez

---

## Overview

This repository contains the code and analysis pipeline for an RNA-sequencing project focused on **ribosome profiling (translatomics)**. The goal is to investigate how the knockout of an unknown gene in human cells affects gene expression at the translational level, as well as ribosomal codon dynamics.

Four samples are analysed: two biological replicates each for wildtype (`WT_Rep1`, `WT_Rep2`) and gene-knockout conditions (`KO_Rep1`, `KO_Rep2`). The analysis covers differential expression, Gene Ontology enrichment, and codon occupancy.

---

## Repository Structure

```
rna-seq/
├── 01_preprocessing/       # Quality control and adapter trimming scripts
├── 02_contamination/       # Removal of undesired RNAs (rRNA, tRNA, snoRNA, snRNA)
├── 03_mapping/             # Genome and transcriptome alignment
├── 04_qc_alignment/        # Alignment quality control (reading frame analysis)
├── 05_read_counting/       # Read count tables per gene and biotype
├── 06_differential_expression/  # DESeq2-based differential expression analysis
├── 07_gene_ontology/       # Gene Ontology enrichment analysis
├── 08_codon_occupancy/     # A-site codon occupancy analysis
└── README.md
```

> **Note:** Adjust folder names above to match the actual structure of your repository.

---

## Analysis Pipeline

### 1. Library Processing & Quality Control

Raw reads were assessed for quality using **FastQC** and summarised with **MultiQC**. Reads were found to contain Illumina Universal Adapters and a four-nucleotide random sequence at the 3′-end characteristic of ribosome profiling. These were removed using **Cutadapt**, and reads shorter than 25 nt were discarded.

| Tool     | Version |
|----------|---------|
| FastQC   | 0.11.9  |
| MultiQC  | 1.8     |
| Cutadapt | 2.5     |

### 2. Removal of Contaminant RNAs

Processed reads were mapped against a combined FASTA file of undesired human RNAs using **Bowtie** to remove contaminants prior to downstream analysis. The contaminant reference was built from:

- sno-RNA, sn-RNA, and rRNA sequences from **Ensembl BioMart**
- tRNA sequences from the **Genomic tRNA Database (GtRNAdb)**
- Additional human rRNAs from the **NCBI Nucleotide Database**

| Tool   | Version |
|--------|---------|
| Bowtie | 1.2.0   |

### 3. Genome Mapping

The filtered read library was mapped to the human reference genome using **Bowtie**. The resulting SAM file was converted to a sorted BAM file using **SAMtools**.

- **Reference genome:** `Homo_sapiens.GRCh38.dna.primary_assembly.fa` (Ensembl)

| Tool     | Version |
|----------|---------|
| Bowtie   | 1.2.0   |
| SAMtools | 1.10    |

### 4. Alignment Quality Control

Reading frame distributions were assessed using the R package **Ribo-seQC** to confirm expected ribosome footprint characteristics (30–31 nt fragments predominantly in the 0-reading frame).

| Tool      | Version |
|-----------|---------|
| Ribo-seQC | 0.99.0  |

### 5. Read Counting

Read count tables were generated using **Subread (featureCounts)**: one for raw gene-level counts and one stratified by RNA biotype (protein-coding, lncRNA, miRNA, etc.).

| Tool    | Version |
|---------|---------|
| Subread | 2.0.1   |

### 6. Differential Expression Analysis

Differential expression between wildtype and knockout was analysed in R using **DESeq2**, with gene annotation provided by **AnnotationDbi** and parallelisation via **BiocParallel**.

| R Package      | Notes                          |
|----------------|--------------------------------|
| DESeq2         | Differential expression        |
| BiocParallel   | Parallelised computation       |
| AnnotationDbi  | Gene annotation                |

### 7. Gene Ontology Analysis

Differentially expressed genes were clustered by biological process, molecular function, and cellular component using **topGO** and **clusterProfiler** in R.

| R Package      | Notes                          |
|----------------|--------------------------------|
| topGO          | GO enrichment (Fisher's test)  |
| clusterProfiler| Enrichment visualisation       |
| AnnotationDbi  | GO term annotation             |

### 8. Codon Occupancy Analysis

Processed reads were additionally mapped to the annotated human transcriptome (extended by 18 nt in both 5′ and 3′ directions to capture footprints at ORF boundaries) using **Bowtie**. A-site codon occupancy was then calculated using the bash script from [LeidelLab/Codon_occupancy_cal](https://github.com/LeidelLab/Codon_occupancy_cal) and visualised in R.

- **Reference transcriptome:** `GRCh38_p13_APPRIS_CDS_plus18` (Ensembl BioMart, ±18 nt extended)

---

## Reference Data

| Data                        | Source                          |
|-----------------------------|---------------------------------|
| Human genome (GRCh38)       | Ensembl                         |
| Human transcriptome (GRCh38_p13_APPRIS_CDS_plus18) | Ensembl BioMart |
| snoRNA / snRNA / rRNA       | Ensembl BioMart                 |
| tRNA                        | Genomic tRNA Database (GtRNAdb) |
| Additional rRNA             | NCBI Nucleotide Database        |

---

## Key Results

- ~63% of reads mapped to the human genome after contamination removal; ~91% of those mapped to protein-coding genes.
- PCA showed greater variance between replicates than between genotypes, indicating a significant and non-random effect of the knockout on gene expression.
- The knockout broadly upregulated stress-related and apoptotic pathways (notably *BMP2*, *PTN*, *FGF13*) and downregulated biosynthetic RNA/ribosomal pathways (notably *MAGEB2*, *UNCX*, *NEFM*).
- Gene Ontology analysis pointed towards synovial fibroblast origin, with patterns consistent with increased environmental sensing and impaired biosynthetic capacity.
- Average codon occupancy decreased by ~3.1% after knockout (median: 8.6% faster), with AGA showing the largest individual shift (+68%, i.e. 40% slower translation). Results were statistically significant (sign-test, p = 8.7E-4).

---

## Dependencies

### Command-line tools
- FastQC ≥ 0.11.9
- MultiQC ≥ 1.8
- Cutadapt ≥ 2.5
- Bowtie ≥ 1.2.0
- SAMtools ≥ 1.10
- Subread ≥ 2.0.1

### R packages
- Ribo-seQC (0.99.0)
- DESeq2
- BiocParallel
- AnnotationDbi
- topGO
- clusterProfiler

### External scripts
- [LeidelLab/Codon_occupancy_cal](https://github.com/LeidelLab/Codon_occupancy_cal)

---

## Usage

Each numbered directory contains scripts for its respective analysis step. Scripts are intended to be run in order. See the comments within each script for parameter details and expected input/output files.

---

## Course Context

This project was completed as part of the RNA-seq module of the MSc Bioinformatics & Computational Biology programme at Universität Bern (2022). The full written report is available separately.
