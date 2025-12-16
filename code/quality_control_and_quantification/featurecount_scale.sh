#!/bin/bash
#本脚本用于对STAR比对生成的BAM文件进行featureCounts计数

BAM_DIR="/public1/home/a8s000332/data_space/RNA_zy/spikein/9_code/05_alignment"
OUT_DIR="/public1/home/a8s000332/data_space/RNA_zy/spikein/9_code/06_STAR/01featurecount/result"
GTF="/public1/home/a8s000332/data_space/RNA_zy/spikein/88_ref/merged_annotation.gtf"

for BAM_FILE in "$BAM_DIR"/*.sorted.bam
do 
  SAMPLE_NAME=$(basename "$BAM_FILE" .sorted.bam)
  echo "Processing:$SAMPLE_NAME"
  featureCounts -p --countReadPairs -T 8 -s 2 -J -B -C -M --fraction\
    -a "$GTF" \
    -o "$OUT_DIR/${SAMPLE_NAME}.counts.txt" \
    "$BAM_FILE"
done
