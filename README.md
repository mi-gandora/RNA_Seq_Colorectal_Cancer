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
* **FeatureCounts**: Alignment and quantification.
* **R (DESeq2/edgeR)**: Statistical analysis and visualization.
