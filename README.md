# Transcriptomic Analysis of Colorectal Cancer (CRC)
## Differential Gene Expression: Tumor vs. Matched Adjacent Normal Tissue

### Project Overview
This project investigates the gene expression profiles of colorectal cancer patients. By comparing tumor tissue with matched adjacent normal tissue from the same individuals, we aim to identify significantly differentially expressed genes (DEGs) while controlling for inter-patient genetic variability.

### Objectives
* Perform Quality Control (QC) on raw RNA-Seq reads.
* Align reads to the human reference genome (GRCh38).
* Quantify gene-level counts.
* Conduct DGEA to identify up-regulated and down-regulated genes in tumor samples.

### Tools Used
* **Linux/Bash**: Data management and preprocessing.
* **FastQC**: Quality assessment.
* **fastp**: Ultra-fast all-in-one FASTQ preprocessing and adapter trimming.
* **BWA (v0.7.17)**: Burrows-Wheeler Aligner for short-read reference mapping.
* **Samtools (v1.x)**: Genomic alignment manipulation, streaming, and coordinate sorting.
* **FeatureCounts**: Post-alignment gene-level quantification.
* **R (DESeq2/edgeR)**: Downstream statistical analysis and visualization.

---

## Pipeline & Step-by-Step Methodology

The complete upstream workflow is automated and documented in the master shell script: `download_samples.sh`. The pipeline processes raw sequencing data sequentially through the following stages:

### 1. Data Acquisition (NCBI SRA)
Raw single-end sequencing data for the tumor and matched normal samples were fetched directly from the NCBI Sequence Read Archive (SRA) database using the SRA Toolkit. 

### 2. Initial Quality Control (FastQC)
Initial raw sequence quality assessment was performed using `FastQC` to evaluate per-base sequence quality, GC content, duplication levels, and detect adapter contamination before any downstream modifications.

### 3. Smart Trimming & Adapter Removal (fastp)
To remove low-quality bases and technical artifacts, `fastp` was utilized with the following parameters:
* Automated adapter detection for sequencing adapters.
* Sliding window trimming at the front and tail ends (`--cut_front`, `--cut_tail`).
* Minimum length filtering (`--length_required 36`) to discard reads compromised by heavy trimming.
* Generation of an interactive HTML quality snapshot (`_fastp_report.html`) for post-trimming validation.

### 4. Memory-Optimized Reference Alignment (BWA-MEM) & Sorting
Cleaned reads were mapped against the human reference genome (**GRCh38/hg38**), followed by coordinate sorting, indexing, and mapping quality assessment.

## Phase 1 Final Artifacts (Alignment & Preprocessing)

The upstream preprocessing pipeline automated by `download_samples.sh` successfully concludes with fully optimized, coordinate-sorted, and indexed BAM files. These serve as the finalized inputs required for the next phase of the project (**Phase 2: Gene Quantification via `featureCounts` and Downstream DGEA**).

```text
├── S1_Tumor_sorted.bam       # Coordinate-aligned BAM file (Input for featureCounts)
├── S1_Tumor_sorted.bam.bai   # Index file for rapid coordinate lookup
│
├── S2_Tumor_sorted.bam       # Coordinate-aligned BAM file (Input for featureCounts)
├── S2_Tumor_sorted.bam.bai   # Index file for rapid coordinate lookup
│
├── S3_Tumor_sorted.bam       # Coordinate-aligned BAM file (Input for featureCounts)
├── S3_Tumor_sorted.bam.bai   # Index file for rapid coordinate lookup
│
├── S4_Tumor_sorted.bam       # Coordinate-aligned BAM file (Input for featureCounts)
└── S4_Tumor_sorted.bam.bai   # Index file for rapid coordinate lookup
│
├── S5_Tumor_sorted.bam       # Coordinate-aligned BAM file (Input for featureCounts)
└── S5_Tumor_sorted.bam.bai   # Index file for rapid coordinate lookup
