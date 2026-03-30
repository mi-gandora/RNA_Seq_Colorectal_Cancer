#!/bin/bash
# 1. Define the Accession ID
ACC="ERR14837534"
echo "-------------------------------------------------------"
echo "Starting Pipeline for: $ACC"
echo "-------------------------------------------------------"

# 2. Secure Download (Handles resume if internet blinks)
echo "[STEP 1/3] Downloading SRA bundle..."
prefetch $AC
C
# 3. Extract FASTQ (Local process, no internet needed)
echo "[STEP 2/3] Extracting FASTQ file..."
fasterq-dump --threads 4 --temp ./tmp_dir $ACC && rm -rf ./tmp_dir

# 4. Immediate Compression (Saves disk space on your 8GB laptop)
echo "Compressing raw file to save space..."
pigz "${ACC}.fastq"

# 5. Initial Quality Control (The 'Before' Snapshot)
# FastQC reads the .gz file directly!
echo "[STEP 4/5] Running FastQC on compressed file..."
fastqc --memory 2048 "${ACC}.fastq.gz"

# 6. Smart Trimming & Verification (The 'After' Snapshot)
# fastp generates its own 'After' report automatically
echo "[STEP 5/5] Trimming and generating Final Report..."
fastp -i "${ACC}.fastq.gz" \
      -o "${ACC}_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_tail \
      --length_required 36 \
      --html "${ACC}_fastp_report.html"

echo "-------------------------------------------------------"
echo "PROCESS COMPLETE"
echo "Raw QC: ${ACC}_fastqc.html"
echo "Trimmed QC: ${ACC}_fastp_report.html"
echo "Clean Data: ${ACC}_trimmed.fastq.gz"
echo "-------------------------------------------------------"
