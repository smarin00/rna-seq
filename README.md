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
├── scripts/
│   ├── 0_download_sra.slurm            # Download raw reads from SRA
│   ├── 1_fast_qc.slurm                 # Pre-processing quality control (FastQC)
│   ├── 2_preprocessing.slurm           # Adapter trimming (Cutadapt)
│   ├── 3_get_indices.slurm             # Build Bowtie indices for contaminant & genome references
│   ├── 4_remove_unwantedRNA.slurm      # Remove contaminant RNAs (rRNA, tRNA, snoRNA, snRNA)
│   ├── 5_fastqc_postprocessing.slurm   # Post-trimming quality control (FastQC)
│   ├── 6_genome_mapping.slurm          # Align reads to human genome (GRCh38)
│   ├── 7_DESeq_counts.slurm            # Generate count tables & run DESeq2 + Gene Ontology
│   ├── 8_transcriptome_mapping.slurm   # Align reads to human transcriptome
│   ├── 9_Codon_Occupancy.sh            # A-site codon occupancy analysis
│   ├── 10_cleanup.slurm                # Remove intermediate files
│   ├── Codon_occupancy_cal.sh          # Codon occupancy calculation helper
│   ├── module.sh                       # Environment module loader
│   ├── module2.sh                      # Environment module loader (alternative)
│   ├── samples.txt                     # Sample IDs list
│   └── samples_name.txt                # Sample display names
└── README.md
```

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

All scripts live in the `scripts/` folder and are numbered to reflect execution order. They are written as SLURM batch scripts (`.slurm`) for HPC cluster submission, with the exception of the codon occupancy steps which are plain bash (`.sh`). Run them sequentially:

```bash
sbatch scripts/0_download_sra.slurm
sbatch scripts/1_fast_qc.slurm
sbatch scripts/2_preprocessing.slurm
# ... and so on through 10_cleanup.slurm
```

Sample IDs and names are defined in `scripts/samples.txt` and `scripts/samples_name.txt`. Module loading for HPC environments is handled by `module.sh` and `module2.sh`.

---

## Course Context

This project was completed as part of the RNA-seq module of the MSc Bioinformatics & Computational Biology programme at Universität Bern (2022). The full written report is available separately.
