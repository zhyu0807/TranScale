import pandas as pd
import glob
import os
#本脚本用于合并多个 featureCounts 生成的 .counts.txt 文件
# 路径设置
count_dir = "/public1/home/a8s000332/data_space/RNA_zy/spikein/9_code/06_STAR/01featurecount/result"
output_file = "/public1/home/a8s000332/data_space/RNA_zy/spikein/9_code/06_STAR/01featurecount/FeatureCounts.gene_quant.csv"

# 获取所有 .counts.txt 文件
count_files = glob.glob(os.path.join(count_dir, "*.counts.txt"))

# 初始化
merged_df = None

for file in count_files:
    # 提取样本名（去掉路径和扩展名）
    sample_name = os.path.basename(file).replace(".counts.txt", "")
    
    # 读取文件，跳过注释（以 # 开头的行）
    df = pd.read_csv(file, sep="\t", comment='#')
    
    # 保留 Geneid 和该样本的 counts 值
    df = df[["Geneid", df.columns[-1]]]
    df.columns = ["Ensembl", sample_name]

    # 合并
    if merged_df is None:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on="Ensembl", how="outer")

# 将 NaN 填为 0（有些基因在某些样本中可能缺失）
merged_df = merged_df.fillna(0).astype({col: int for col in merged_df.columns if col != "Ensembl"})

# 保存为 CSV
merged_df.to_csv(output_file, index=False)
print(f" 表达矩阵已保存至：{output_file}")

