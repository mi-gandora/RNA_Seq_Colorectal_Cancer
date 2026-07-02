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
* **R (edgeR)**: Downstream statistical analysis and visualization.

---
# Phase 1: Upstream Processing (Linux / WSL)

## Pipeline & Step-by-Step Methodology

The complete upstream workflow is automated and documented in the master shell script: `download_samples.sh`. The pipeline processes raw sequencing data sequentially through the following stages, concluding with optimized gene-level quantification:

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

### 5. Gene Quantification & Pipeline Optimization
To transform the raw sequence alignments into a structured gene expression profile, all 10 coordinate-sorted BAM files (5 Tumor + 5 Normal pairs) were quantified simultaneously using `featureCounts`.

### Phase 1 Directory Output

The terminal pipeline successfully concludes here. While the large intermediate alignment files remain archived locally in Linux as raw biological evidence, the entire upstream workflow has been compressed into a single, production-ready spreadsheet for statistical analysis.

```text
├── Intermediate Data (Archived locally in Linux)
│   ├── *_sorted.bam          # 10 Coordinate-aligned sequence files (5 Tumor / 5 Normal)
│   └── *_sorted.bam.bai      # 10 Genomic index maps for rapid coordinate lookup
│
└── Final Production Output (The entry point for Phase 2)
    └── counts.csv            # Consolidated Master Count Matrix (Sole input for R)
```

# Phase 2: Downstream Statistical Analysis & Functional Interpretation (R)

## Pipeline & Step-by-Step Methodology
The downstream analysis is executed using specialized R scripts designed to parse the quantified master count matrix `(counts.csv)`, isolate confident statistical variations, and map them to known physiological pathways.

### 1. Paired Experimental Design Setup (edgeR/limma)
To effectively isolate cancer-driving signals from natural human genetic variation, a paired-sample linear model was implemented. By tracking both the clinical condition (Tumor vs. Adjacent Normal) and individual patient origins (P1 through P5), the model subtracts inter-patient baseline differences to maximize mathematical power.

### 2. Differential Gene Expression Analysis (DGEA)
* Statistical Threshold: Adjusted p-value calculated using the Benjamini-Hochberg (BH) False Discovery Rate (FDR) multiple-testing correction.  
* Biological Threshold: A minimum 2-fold expression difference (log2FC > 1.0) to isolate high-confidence targets. 
* Annotation & Verification: Ensembl IDs were mapped to official HGNC symbols using biomaRt. Missing gene symbols default to their stable Ensembl IDs, and biological duplicates are safely resolved using `make.unique()` to preserve downstream matrix structure.

### 3. Transcriptional Profile Visualization
A global Volcano Plot was generated to map the biological effect size against statistical confidence. By adjusting for the matched patient backgrounds, background noise is thoroughly suppressed, causing highly confident transcripts to stand out clearly
* Upregulated: Significant genes overexpressed in the colorectal tumor matrix (log2FC > 1.0, p_adj < 0.05).
* Downregulated: Significant genes suppressed in the colorectal tumor matrix (log2FC < -1.0, p_adj < 0.05).
* Insignificant: Transcripts failing to clear the combined biological or statistical checkpoints

### 4. Over-Representation Analysis (ORA) via Gene Ontology
Significant DEGs were mapped to their functional coordinates across Biological Process (BP), Cellular Component (CC), and Molecular Function (MF) domains using `enrichGO()`. Enriched terms were programmatically streamlined using `clusterProfiler::simplify()` to eliminate redundant parental vocabulary and isolate the core active biological themes.

### 5. Gene Set Enrichment Analysis (GSEA)
To capture coordinated, system-wide pathway shifts that standard strict cutoff methods might obscure, a threshold-free GSEA model was executed via `gseGO()`. Every single quantified gene was retained and continuously ranked based on its log-fold change expression gradient (log2FC).

### Phase 2 Directory Output
The complete downstream analysis scripts, curated spreadsheets, and visual plots generated in Phase 2 are structured within the repository as follows:

```text
├── Scripts
│   ├── Differential Gene Expression Analysis.R  # Paired limma-voom script, annotations, and GO profiles
│   └── Gene Set Enrichment Analysis_2.R       # Threshold-free continuous GSEA pipeline script
│
├── Data_Outputs
│   ├── total_degs.csv                          # Polished spreadsheet containing all identified DEGs
│   ├── upregulated_genes.csv                   # Subset of features significantly turned ON in CRC
│   ├── downregulated_genes.csv                 # Subset of features significantly turned OFF in CRC
│   └── GSEA_result.csv                         # Full statistical table of enriched gene set terms
│
└── Visualizations
    ├── volcanoplots.png                        # Differential Expression landscape plot
    ├── Bp dotplot.png / CC dotplot.png         # Simplified GO ontology term distributions
    ├── gse_dotplot.png                         # GSEA dotplot showing activated and suppressed pathways
    └── gse1.png / gse2.png                     # Enrichment running plots for target cascades
=======



