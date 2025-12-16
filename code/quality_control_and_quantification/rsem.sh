#!/bin/bash

FASTQ_DIR="/public1/home/a8s000332/data_space/RNA_zy/spikein/01_clean_fq"
REF_PREFIX="/public1/home/a8s000332/data_space/RNA_zy/spikein/rsem/ref/human_merged"
OUT_DIR="/public1/home/a8s000332/data_space/RNA_zy/spikein/rsem/result"
THREADS=16

mkdir -p "$OUT_DIR"

# 链特异性
STRANDED="reverse"    

for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    echo "===== Processing $SAMPLE ====="

    rsem-calculate-expression \
        --paired-end \
        --star \
	--star-gzipped-read-file \
        --strandedness $STRANDED \
        --star-output-genome-bam \
        --estimate-rspd \
        --append-names \
        "$R1" "$R2" \
        "$REF_PREFIX" \
        "${OUT_DIR}/${SAMPLE}" \
        --num-threads ${THREADS}

    echo "===== Done $SAMPLE ====="
done

