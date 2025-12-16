# -*- coding: utf-8 -*-
#检索指定文件夹下面的fq.gz文件并生成一个sample Sample_name为表头的文件
#在每次运行前需要确定样本后缀是否为 _1.fq.gz和 _2.fq.gz，不是的话需要对代码进行修改
import os 
import pandas as pd 
import argparse

def generate_fq_list(input_dir,output_file):
    
    if not os.path.exists(input_dir):
        print("The inut dir {input_dir} is not exist")
   
    sample_info = {}
   
    with open(output_file,'w') as f:
        f.write("Sample\tSample_path_q1\tSample_path_q2\n")
    
    for root,dirs,files in os.walk(input_dir):
        for file in files:
            
            if file.endswith('.fq.gz'):
                absolute_path = os.path.abspath(os.path.join(root,file))
                
                if file.endswith('_1.fq.gz'):
                    Sample_name = file.split('_1.fq.gz')[0]
                    
                    if Sample_name not in sample_info:
                        sample_info[Sample_name] = {"Sample_path_q1":"", "Sample_path_q2":""}
                    sample_info[Sample_name]["Sample_path_q1"] = absolute_path
                elif file.endswith('_2.fq.gz'):
                    Sample_name = file.split('_2.fq.gz')[0]
                    if Sample_name not in sample_info:
                        sample_info[Sample_name] = {"Sample_path_q1":"", "Sample_path_q2":""}
                    sample_info[Sample_name]["Sample_path_q2"] = absolute_path
    
    df = pd.DataFrame(sample_info).T
    df.reset_index(inplace=True)
    df.rename(columns={'index':'Sample'},inplace=True)
    
    df.to_csv(output_file,sep='\t',index=False)
    
    print(f"文件列表已生成，文件路径为{output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="生成fq.gz文件列表")
    parser.add_argument("-i","--input_dir",help="输入文件夹路径",required=True)
    parser.add_argument("-o","--output_file",help="输出文件路径",required=True)
    args = parser.parse_args()
    generate_fq_list(args.input_dir,args.output_file)
