import pandas as pd
import glob
import os
import re

#本脚本为了合并多个 RSEM 生成的 .genes.results 文件，最终生成excle格式的表格文件
files = glob.glob("/public1/home/a8s000332/data_space/RNA_zy/spikein/rsem/result/*.genes.results")
files.sort() 

if not files:
    print("错误：未找到 .results 文件")
    exit()

print(f"找到 {len(files)} 个文件，开始处理...")


df_counts = pd.DataFrame()
df_tpm = pd.DataFrame()
df_fpkm = pd.DataFrame()

for file in files:
    
    sample_name = file.split(".sorted_spikein")[0]
    
    
    data = pd.read_csv(file, sep="\t", index_col="gene_id")
    
    
    if df_counts.empty:
        df_counts = pd.DataFrame(data['expected_count']).rename(columns={'expected_count': sample_name})
        df_tpm = pd.DataFrame(data['TPM']).rename(columns={'TPM': sample_name})
        df_fpkm = pd.DataFrame(data['FPKM']).rename(columns={'FPKM': sample_name})
    else:
        df_counts[sample_name] = data['expected_count']
        df_tpm[sample_name] = data['TPM']
        df_fpkm[sample_name] = data['FPKM']


df_counts.reset_index(inplace=True)
df_tpm.reset_index(inplace=True)
df_fpkm.reset_index(inplace=True)


print("正在写入 Excel...")
df_counts.to_excel("STAR_Hista_Counts.xlsx", index=False)
df_tpm.to_excel("STAR_Hista_TPM.xlsx", index=False)
df_fpkm.to_excel("STAR_Hista_FPKM.xlsx", index=False)

print("完成！文件已保存在当前目录。")
