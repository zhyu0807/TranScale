import os
import pandas as pd

# 加载样本及样本信息
configfile: "bam_config.yaml"

SAMPLES = pd.read_csv(config["samples_file"], sep="\t")
SAMPLES = SAMPLES.rename(columns={
    "Sample": "sample",
    "Sample_path_q1": "q1",
    "Sample_path_q2": "q2"
}).set_index("sample", drop=False)

rule all:
    input:
        expand("03_quality_control/qualimap/{sample}_report/qualimapReport.html", sample=SAMPLES.index),
        expand("05_alignment/{sample}.sorted.bam", sample=SAMPLES.index),
        expand("05_alignment/{sample}.sorted.bam.bai", sample=SAMPLES.index)
rule star_align:
    input:
        r1 = lambda wc: SAMPLES.loc[wc.sample, "q1"],
        r2 = lambda wc: SAMPLES.loc[wc.sample, "q2"],
        index = config["reference"]["star_index_dir"]
    output:
        temp("05_alignment/{sample}.Aligned.sortedByCoord.out.bam"),
        "05_alignment/{sample}.ReadsPerGene.out.tab"
    params:
        prefix = lambda wc: f"05_alignment/{wc.sample}.",
        tmp_dir = lambda wc: f"/public1/home/a8s000332/data_space/RNA_zy/RNA_spike_in/tmp/{wc.sample}_STARtmp",
        threads = lambda wc: config["threads"]["star_align"]
    resources:
        star = 1
    log:
        "logs/star_align/{sample}.log"
    shell:
        """
        STAR --runThreadN {params.threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outTmpDir {params.tmp_dir} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --twopassMode Basic \
             --quantMode GeneCounts \
             --outSAMunmapped Within > {log} 2>&1
        """

rule index_bam:
    input:
        "05_alignment/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam = "05_alignment/{sample}.sorted.bam",
        bai = "05_alignment/{sample}.sorted.bam.bai"
    log:
        "logs/index_bam/{sample}.log"
    shell:
        """
        # 将STAR输出的默认文件名重命名为更简洁的最终文件名
        mv {input} {output.bam}
        # 为重命名后的BAM文件创建索引，samtools会自动生成.bai文件
        samtools index {output.bam} > {log} 2>&1
        """

rule qualimap:
    input:
        bam = "05_alignment/{sample}.sorted.bam",
        bai = "05_alignment/{sample}.sorted.bam.bai",
        gtf = config["reference"]["annotation"]
    output:
        report = "03_quality_control/qualimap/{sample}_report/qualimapReport.html"
    params:
        mem = lambda wc: config["qualimap"]["java_mem"],
        output_dir = lambda wc: f"03_quality_control/qualimap/{wc.sample}_report"
    log:
        "logs/qualimap/{sample}.log"
    shell:
        """
        export _JAVA_OPTIONS="-Djava.awt.headless=true -Xmx{params.mem}"
        mkdir -p {params.output_dir}
        qualimap rnaseq -bam {input.bam} \
                       -gtf {input.gtf} \
                       -outformat PDF:HTML \
                       -outdir {params.output_dir} > {log} 2>&1
        """
