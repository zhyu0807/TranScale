import os
import pandas as pd

# 设置参数
work_path = "/public1/home/a8s000332/data_space/RNA_zy/spikein"

# 读取样本信息
SAMPLES = pd.read_csv(f"{work_path}/9_code/sample.txt", sep="\t")
SAMPLES.rename(columns={
    "Sample": "sample",
    "Sample_path_q1": "q1",
    "Sample_path_q2": "q2"
}, inplace=True)

INPUT_PATH = f"{work_path}/00_raw_fq"
CLEAN_DATA = f"{work_path}/01_clean_fq"
FASTQC_REPORTS = f"{work_path}/01_fastqc_reports"

# 创建目录
if not os.path.exists(CLEAN_DATA):
    os.makedirs(CLEAN_DATA)
if not os.path.exists(FASTQC_REPORTS):
    os.makedirs(FASTQC_REPORTS)

sample_names = SAMPLES["sample"].tolist()

def get_fastq(sample, read_col):
    return SAMPLES.loc[SAMPLES["sample"] == sample, read_col].values[0]

rule all:
    input:
        expand(CLEAN_DATA + "/{sample}_R1.fastq.gz", sample=sample_names),
        expand(CLEAN_DATA + "/{sample}_R2.fastq.gz", sample=sample_names),
        expand(FASTQC_REPORTS + "/{sample}.tar.gz", sample=sample_names)

rule fastp:
    input:
        r1 = lambda wc: get_fastq(wc.sample, "q1"),
        r2 = lambda wc: get_fastq(wc.sample, "q2")
    output:
        cleaned_r1 = CLEAN_DATA + "/{sample}_R1.fastq.gz",
        cleaned_r2 = CLEAN_DATA + "/{sample}_R2.fastq.gz",
        html = work_path + "/01_fastp_reports/{sample}_fastp.html",
        json = work_path + "/01_fastp_reports/{sample}_fastp.json"
    shell:
        """
        fastp \
        -i {input.r1} \
        -I {input.r2} \
        -o {output.cleaned_r1} \
        -O {output.cleaned_r2} \
        -h {output.html} \
        -j {output.json}
        """

rule fastqc:
    input:
        r1 = CLEAN_DATA + "/{sample}_R1.fastq.gz",
        r2 = CLEAN_DATA + "/{sample}_R2.fastq.gz"
    output:
        out_dir = directory(FASTQC_REPORTS + "/{sample}"),
        out_file = FASTQC_REPORTS + "/{sample}.tar.gz"
    shell:
        """
        mkdir -p {output.out_dir}
        fastqc -o {output.out_dir} {input.r1} {input.r2}
        tar -zcvf {output.out_file} -C {output.out_dir} .
        """
