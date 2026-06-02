#!/bin/bash
# 1. Define the Accession ID
ACC="ERR14837534"
echo "-------------------------------------------------------"
echo "Starting Pipeline for: $ACC"
echo "-------------------------------------------------------"

# 2. Secure Download
echo "[STEP 1/3] Downloading SRA bundle..."
prefetch $AC
C
# 3. Extract FASTQ
echo "[STEP 2/3] Extracting FASTQ file..."
fasterq-dump --threads 4 --temp ./tmp_dir $ACC && rm -rf ./tmp_dir

# 4. Immediate Compression
echo "Compressing raw file to save space..."
pigz "${ACC}.fastq"

# 5. Initial Quality Control
# FastQC reads the .gz file directly
echo "[STEP 4/5] Running FastQC on compressed file..."
fastqc --memory 2048 "${ACC}.fastq.gz"

# 6. Smart Trimming & Verification
echo "[STEP 5/5] Trimming and generating Final Report..."
fastp -i "${ACC}.fastq.gz" \
      -o "${ACC}_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_tail \
      --length_required 36 \
      --html "${ACC}_fastp_report.html"

# 7. Optimized Reference Alignment & Sorting (Resource-Constrained Configuration)
# Environment Restrictions: 8GB RAM / Strict Windows Host Disk Limits
echo "[STEP 6] Running memory-optimized alignment and sorting..."

# Ensure the host temporary directory exists
mkdir -p /mnt/c/bio_tmp

# Low-RAM Streaming Pipeline:
# - bwa mem maps reads against the GRCh38 human reference genome.
# - samtools view -u streams uncompressed BAM to prevent a massive 40GB+ SAM file layout.
# - samtools sort strictly caps memory at 500M to prevent kernel OOM (Out Of Memory) panics.
# - -T redirects cache chunks to the external C: drive to protect the Linux VHDX from inflating.
bwa mem -t 2 ../../reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa "${ACC}_trimmed.fastq.gz" | \
samtools view -u - | \
samtools sort -@ 1 -m 500M -T /mnt/c/bio_tmp/${ACC}_sort -o "${ACC}_sorted.bam" -

# Clean up host disk space immediately after sample success
rm -f /mnt/c/bio_tmp/${ACC}_sort.*.bam

# Index the sorted BAM file
echo "Generating genomic index (.bai) for ${ACC}_sorted.bam..."
samtools index "${ACC}_sorted.bam"

# Alignment Quality Control (Flagstat)
echo "Generating mapping statistics for ${ACC}_sorted.bam..."
samtools flagstat "${ACC}_sorted.bam" > "${ACC}_alignment_stats.txt"

echo "-------------------------------------------------------"
echo "PROCESS COMPLETE"
echo "Raw QC: ${ACC}_fastqc.html"
echo "Trimmed QC: ${ACC}_fastp_report.html"
echo "Clean Data: ${ACC}_trimmed.fastq.gz"
echo "Aligned & Sorted Reads: ${ACC}_sorted.bam"
echo "Genomic Index File: ${ACC}_sorted.bam.bai"
echo "Alignment Quality Report: ${ACC}_alignment_stats.txt"
echo "-------------------------------------------------------"


# ==============================================================================
# PHASE 1 FINAL STAGE: MULTI-SAMPLE GENE QUANTIFICATION
# ==============================================================================
echo "======================================================="
echo "[STEP 7] RUNNING OPTIMIZED GENE QUANTIFICATION"
echo "======================================================="

GTF_REF="../../reference/Homo_sapiens.GRCh38.115.gtf"

# Run featureCounts simultaneously across all 10 generated BAM files
# -s 0        : Unstranded library kit (enforcing -s 1 or -s 2 drops assignment to ~6.7%)
# -t gene     : Gene locus resolution captures overlapping unstranded transcript fragments
# --maxMOp 100 : Expands CIGAR operations to retain highly spliced fragments
featureCounts -T 2 -s 0 -t gene -g gene_id --maxMOp 100 \
  -a "$GTF_REF" \
  -o ./counts_raw.txt \
  *_sorted.bam

echo "-------------------------------------------------------"
echo "[STEP 8] CONVERTING RAW OUTPUT TO PRODUCTION FORMAT"
echo "-------------------------------------------------------"

# Clean raw text matrix into a true comma-separated format for downstream R ingestion
grep -v '^#' ./counts_raw.txt | tr '\t' ',' > ./counts.csv

# Delete intermediate text summaries to maintain a pristine workspace
rm ./counts_raw.txt

echo "=============================================================================="
echo "PHASE 1 COMPLETE: All samples processed. Master matrix 'counts.csv' generated."
echo "Ready for Phase 2 Downstream Analysis in R Studio."
echo "=============================================================================="
