#!/usr/bin/env bash
# =============================================================================
# 01_qc_cutadapt.sh  (shell, not R: these steps run outside R)
# Read QC (FastQC v0.12.1, MultiQC v1.14) and primer removal (Cutadapt v4.4).
#
# Expects raw paired FASTQ from SRA in data/raw/<DATASET>/<RUN>_1.fastq.gz and
# <RUN>_2.fastq.gz (the naming used by fasterq-dump --split-files).
# Writes data/trimmed/<DATASET>/<RUN>_R1.trim.fastq.gz / _R2.trim.fastq.gz
#
# Download example (SRA Toolkit):
#   prefetch SRR... && fasterq-dump --split-files SRR... && gzip SRR*_?.fastq
# =============================================================================
set -euo pipefail

THREADS=8
FWD="CCTACGGGNGGCWGCAG"        # 341F
REV="GACTACHVGGGTATCTAATCC"    # 806R
DATASETS=("SRP411931" "PRJNA822613")

mkdir -p results/qc/fastqc_raw results/qc/cutadapt_logs

# ---- 1. FastQC + MultiQC on raw reads ----------------------------------------
fastqc -t "$THREADS" -o results/qc/fastqc_raw data/raw/*/*.fastq.gz
multiqc -f -o results/qc results/qc/fastqc_raw

# ---- 2. Primer removal ----------------------------------------------------
# --discard-untrimmed : drop read pairs where the primer is not found
# --minimum-length 200: drop reads shorter than 200 bp after trimming
for DS in "${DATASETS[@]}"; do
  mkdir -p "data/trimmed/$DS"
  for R1 in data/raw/"$DS"/*_1.fastq.gz; do
    RUN=$(basename "$R1" _1.fastq.gz)
    R2="data/raw/$DS/${RUN}_2.fastq.gz"
    cutadapt -j "$THREADS" \
      -g "$FWD" -G "$REV" \
      --discard-untrimmed --minimum-length 200 \
      -o "data/trimmed/$DS/${RUN}_R1.trim.fastq.gz" \
      -p "data/trimmed/$DS/${RUN}_R2.trim.fastq.gz" \
      "$R1" "$R2" > "results/qc/cutadapt_logs/${RUN}.log"
  done
done

# ---- 3. Trimming retention summary (the "99.7 %" figure in Section 2.2) ------
echo -e "run\tpairs_in\tpairs_written" > results/qc/cutadapt_summary.tsv
for LOG in results/qc/cutadapt_logs/*.log; do
  RUN=$(basename "$LOG" .log)
  IN=$(grep "Total read pairs processed" "$LOG" | awk '{gsub(",","",$NF); print $NF}')
  OUT=$(grep "Pairs written (passing filters)" "$LOG" | awk '{gsub(",","",$5); print $5}')
  echo -e "$RUN\t$IN\t$OUT" >> results/qc/cutadapt_summary.tsv
done
awk -F'\t' 'NR>1{i+=$2; o+=$3} END{printf "Overall pairs retained after Cutadapt: %.1f%%\n", 100*o/i}' \
  results/qc/cutadapt_summary.tsv
