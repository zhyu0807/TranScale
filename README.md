#TranScale
This repository is used to store the code for data processing and visualization in the article.

Quality Control and Quantification Module
This directory contains the pipeline for initial data processing, quality control, genomic alignment, and transcript quantification.

Workflow Overview
The workflow starts with raw FastQ files, performs quality filtering via fastp, aligns reads using STAR, and generates expression matrices through both featureCounts and RSEM.

File Descriptions
1. Preprocessing & QC
make_filepath.py: A utility to scan directories and generate a sample.txt metadata file with absolute paths to paired-end FastQ files.
RNA_qc.smk: A Snakemake workflow that automates adapter trimming and quality filtering using fastp, followed by sequencing quality assessment via FastQC.

2. Alignment & Mapping QC
RNA_bam.smk: A Snakemake workflow for genomic alignment using the STAR aligner. It includes automated indexing with samtools and post-alignment quality reporting via Qualimap.
bam_config.yaml: The central configuration file defining reference genome paths, tool parameters (e.g., STAR overhang), and resource allocations (CPU/Memory).

3. Quantification
featurecount_scale.sh: A shell script for gene-level quantification using featureCounts, configured for stranded paired-end data and multi-mapping read handling.
featur_count.py: Aggregates individual .counts.txt files from featureCounts into a consolidated "Ensembl-by-Sample" CSV matrix.
rsem.sh: Executes the RSEM pipeline (incorporating STAR) to provide EM-optimized gene and isoform abundance estimates.
result.py: Consolidates RSEM outputs into standardized Excel workbooks containing Counts, TPM, and FPKM metrics for all samples.

If you need the data involved in the project, please contact [Yu Zhang](zhangyu0807@126.com)
