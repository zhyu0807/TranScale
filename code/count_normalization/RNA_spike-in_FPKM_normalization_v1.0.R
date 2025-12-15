# 安装并加载必要的包
if (!require("dplyr")) install.packages("dplyr")
if (!require("edgeR")) install.packages("edgeR")
if (!require("GenomicFeatures")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
}
if (!require("txdbmaker")) {
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install("txdbmaker")
}
if (!require("openxlsx")) {
  install.packages("openxlsx")
}

# 加载包
library(dplyr)
library(edgeR)
library(GenomicFeatures)
library(txdbmaker)
library(openxlsx)

# ========== 第一部分：从count数据计算FPKM ==========

cat("开始第一部分：从count数据计算FPKM\n")

# 1. 读取count矩阵 - 确保所有列都是数值型，按实际替换文件名
count_data <- read.csv("count_data matrix.csv", row.names = 1, header = TRUE, check.names = FALSE)

# 检查数据类型
cat("Count矩阵数据类型:\n")
print(sapply(count_data, class))

# 如果有非数值列，转换为数值型
for (col in colnames(count_data)) {
  if (!is.numeric(count_data[[col]])) {
    cat("将列", col, "转换为数值型\n")
    count_data[[col]] <- as.numeric(as.character(count_data[[col]]))
  }
}

# 查看count矩阵的前几行和维度
cat("Count矩阵维度:", dim(count_data), "\n")
print(head(count_data))

# 2. 使用GenomicFeatures包处理GTF文件并计算外显子并集长度
gtf_file <- "merged_annotation.gtf.gz"

cat("正在使用GenomicFeatures处理GTF文件，这可能需要一些时间...\n")

# 创建TxDb对象
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# 获取每个基因的外显子并集长度
exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths_list <- sum(width(reduce(exons_by_gene)))

# 将列表转换为数据框
gene_lengths_df <- data.frame(
  gene_id = names(gene_lengths_list),
  length = as.numeric(gene_lengths_list),
  stringsAsFactors = FALSE
)

# 输出基因长度列表到CSV文件
write.csv(gene_lengths_df, "gene_lengths.csv", row.names = FALSE)
cat("基因长度列表已保存为 'gene_lengths.csv'\n")

# 查看基因长度数据框
cat("基因长度数据框的前几行（使用外显子并集长度）:\n")
print(head(gene_lengths_df))

# 3. 确保基因ID格式一致并匹配
# 查看count矩阵和基因长度数据的基因ID示例
cat("Count矩阵中的基因ID示例:", head(rownames(count_data), 3), "\n")
cat("GTF文件中的基因ID示例:", head(gene_lengths_df$gene_id, 3), "\n")

# 如果基因ID格式不一致，可能需要进行转换
# 例如，如果count矩阵使用ENSEMBL ID而GTF使用带版本号的ENSEMBL ID
# rownames(count_data) <- gsub("\\..*", "", rownames(count_data))

# 4. 匹配基因长度与count矩阵
# 确保count矩阵中的基因在基因长度数据中都有对应
matched_indices <- match(rownames(count_data), gene_lengths_df$gene_id)
missing_genes <- sum(is.na(matched_indices))

cat("找不到长度信息的基因数量:", missing_genes, "\n")

if (missing_genes > 0) {
  # 移除没有长度信息的基因
  count_data <- count_data[!is.na(matched_indices), ]
  matched_indices <- matched_indices[!is.na(matched_indices)]
  cat("移除缺失基因后，Count矩阵维度:", dim(count_data), "\n")
}

# 提取匹配的基因长度
matched_gene_lengths <- gene_lengths_df$length[matched_indices]

# 5. 计算并验证测序深度处理
# 计算每个样本的总read count（测序深度）
# 确保count_data是数值矩阵
if (!all(sapply(count_data, is.numeric))) {
  cat("警告: count矩阵中仍有非数值列，尝试强制转换\n")
  count_data <- as.data.frame(sapply(count_data, as.numeric))
  rownames(count_data) <- rownames(count_data)  # 保持行名
}

# 清洗 count 矩阵 —— 去掉 NA、负值、非整数
cat("检查并清洗 count 矩阵中的异常值...\n")

# 转成真正的数值矩阵（data.frame 可能带因子）
count_data <- as.matrix(count_data)
mode(count_data) <- "numeric"

# 把所有 NA 或负值设为 0
count_data[is.na(count_data) | count_data < 0] <- 0

# 检查并处理非整数
if (!all(count_data == as.integer(count_data))) {
  cat("检测到非整数 count，四舍五入取整...\n")
  count_data <- round(count_data)
}

# 再检查一次
cat("剩余 NA 数量:", sum(is.na(count_data)), "\n")
cat("Count 矩阵范围:", range(count_data), "\n")

# 6. 计算FPKM (Fragments Per Kilobase Million)
cat("正在计算FPKM...\n")

# 使用edgeR的rpkm函数计算FPKM
dge <- DGEList(counts = count_data)
fpkm_matrix <- rpkm(dge, gene.length = matched_gene_lengths, log = FALSE)

# 查看FPKM矩阵
cat("FPKM矩阵维度:", dim(fpkm_matrix), "\n")
print(head(fpkm_matrix))

# 7. 保存结果
# 将FPKM矩阵转换为数据框，并将基因ID作为第一列
fpkm_df <- as.data.frame(fpkm_matrix)
fpkm_df <- cbind(gene_id = rownames(fpkm_df), fpkm_df)
rownames(fpkm_df) <- NULL

# 保存FPKM矩阵，A1单元格将显示"gene_id"
write.csv(fpkm_df, "fpkm_matrix.csv", row.names = FALSE)
cat("FPKM矩阵已保存为 'fpkm_matrix.csv'，A1单元格为'gene_id'\n")

# 8. 可选：对FPKM值进行log2转换（常用于下游分析）
log2_fpkm_matrix <- log2(fpkm_matrix + 1)
log2_fpkm_df <- as.data.frame(log2_fpkm_matrix)
log2_fpkm_df <- cbind(gene_id = rownames(log2_fpkm_df), log2_fpkm_df)
rownames(log2_fpkm_df) <- NULL
write.csv(log2_fpkm_df, "log2_fpkm_matrix.csv", row.names = FALSE)
cat("log2转换后的FPKM矩阵已保存为 'log2_fpkm_matrix.csv'\n")

cat("\n第一部分处理完成！FPKM计算已完成。\n")

# ========== 第二部分：计算输入拷贝数和拆分及筛选转换fpkm矩阵 ==========

cat("\n\n开始第二部分：计算输入拷贝数和拆分及筛选转换fpkm矩阵\n")

# 读取文件
ref_data <- read.csv("reference value.csv", header = TRUE)
fpkm_matrix <- read.csv("fpkm_matrix.csv", header = TRUE)

# 获取fpkm_matrix文件的数据列数（从第2列开始）
data_cols_count <- ncol(fpkm_matrix) - 1

# 步骤1: 生成input reference value.csv
# 创建新的数据框，第一列与reference文件相同
new_ref_data <- data.frame(ref_data[, 1])
colnames(new_ref_data)[1] <- colnames(ref_data)[1]

# 获取reference文件的B列和C列（假设B列是第2列，C列是第3列）
ref_B <- ref_data[, 2]  # B列
ref_C <- ref_data[, 3]  # C列

# 生成新的数据列
for (i in 1:ceiling(data_cols_count / 3)) {
  # 确定使用B列还是C列（奇数组用B列，偶数组用C列）
  if (i %% 2 == 1) {
    # 奇数组：使用B列
    ref_values <- ref_B
    ref_col_name <- "B"
  } else {
    # 偶数组：使用C列
    ref_values <- ref_C
    ref_col_name <- "C"
  }
  
  # 计算当前组应该生成几列（最后一组可能不足3列）
  start_col <- (i - 1) * 3 + 2
  end_col <- min(start_col + 2, ncol(fpkm_matrix))
  cols_in_this_group <- end_col - start_col + 1
  
  # 生成当前组的列
  for (j in 1:cols_in_this_group) {
    col_index <- start_col + j - 1
    col_name <- colnames(fpkm_matrix)[col_index]
    new_ref_data[[col_name]] <- ref_values * 0.4
  }
  
  cat(sprintf("第%d组: 使用reference %s列，生成%d列 (从%d到%d)\n", 
              i, ref_col_name, cols_in_this_group, start_col, end_col))
}

# 保存input reference value.csv文件
write.csv(new_ref_data, "input reference value.csv", row.names = FALSE)

# 步骤2: 生成log input reference value.csv
# 创建对数转换版本
log_ref_data <- new_ref_data

# 对数据列进行以2为底的对数转换（跳过第一列基因名）
for (col in 2:ncol(log_ref_data)) {
  # 检查是否有非正值，这些值不能取对数
  if (any(log_ref_data[[col]] <= 0)) {
    warning(paste("列", colnames(log_ref_data)[col], "包含非正值，将在对数转换前处理"))
    # 将非正值替换为一个很小的正数，避免对数计算错误
    log_ref_data[[col]][log_ref_data[[col]] <= 0] <- 1e-10
  }
  log_ref_data[[col]] <- log2(log_ref_data[[col]])
}

# 保存log input reference value.csv文件
write.csv(log_ref_data, "log input reference value.csv", row.names = FALSE)

# 步骤3: 从fpkm_matrix.csv提取observed spike-in value.csv
# 按照input reference value.csv数据表中的第一列基因顺序提取

# 获取input reference value.csv的第一列基因名
gene_col_name <- colnames(new_ref_data)[1]
input_ref_genes <- new_ref_data[, 1]

# 确保fpkm_matrix文件中有相同的第一列名称
if (!gene_col_name %in% colnames(fpkm_matrix)) {
  stop("在fpkm_matrix.csv中找不到与input reference value.csv相同的第一列名称")
}

# 初始化observed_spikein_data数据框，使用与input reference value.csv相同的列结构
observed_spikein_data <- data.frame(matrix(nrow = 0, ncol = ncol(new_ref_data)))
colnames(observed_spikein_data) <- colnames(new_ref_data)

# 记录哪些基因被提取了
extracted_genes <- c()

# 提取fpkm_matrix文件中对应的数据
for (i in 1:length(input_ref_genes)) {
  gene <- input_ref_genes[i]
  
  # 在fpkm_matrix文件中找到对应的行
  row_index <- which(fpkm_matrix[, gene_col_name] == gene)
  
  if (length(row_index) == 0) {
    warning(paste("基因", gene, "在fpkm_matrix.csv中未找到"))
    # 如果找不到基因，创建一行NA值
    new_row <- rep(NA, ncol(new_ref_data))
    new_row[1] <- gene  # 第一列是基因名
  } else {
    # 记录被提取的基因
    extracted_genes <- c(extracted_genes, row_index)
    
    # 创建新行，第一列是基因名
    new_row <- c(gene)
    
    # 提取fpkm_matrix文件中对应基因的数据列
    # 确保只提取与input reference value.csv相同列名的数据
    for (col_name in colnames(new_ref_data)[-1]) {
      if (col_name %in% colnames(fpkm_matrix)) {
        new_row <- c(new_row, fpkm_matrix[row_index, col_name])
      } else {
        warning(paste("列", col_name, "在fpkm_matrix.csv中不存在"))
        new_row <- c(new_row, NA)
      }
    }
  }
  
  # 将新行添加到observed_spikein_data中
  observed_spikein_data <- rbind(observed_spikein_data, new_row)
}

# 确保列名正确
colnames(observed_spikein_data) <- colnames(new_ref_data)

# 确保数据类型正确（数据列应该是数值型）
for (col in 2:ncol(observed_spikein_data)) {
  observed_spikein_data[, col] <- as.numeric(observed_spikein_data[, col])
}

# 保存observed spike-in value.csv文件
write.csv(observed_spikein_data, "observed spike-in value.csv", row.names = FALSE)

# 步骤3.5: 创建log observed spike-in value.csv（对数转换的观察值数据）
# 创建对数转换版本
log_observed_spikein_data <- observed_spikein_data

# 对数据列进行以2为底的对数转换（跳过第一列基因名）
if (nrow(log_observed_spikein_data) > 0) {
  for (col in 2:ncol(log_observed_spikein_data)) {
    # 检查是否有0值，如果有0值则加上0.01
    zero_indices <- which(log_observed_spikein_data[[col]] == 0)
    if (length(zero_indices) > 0) {
      cat(sprintf("在列 %s 中发现 %d 个0值，将加上0.01\n", 
                  colnames(log_observed_spikein_data)[col], length(zero_indices)))
      log_observed_spikein_data[[col]][zero_indices] <- 0.01
    }
    
    # 检查是否有负值，这些值不能取对数
    if (any(log_observed_spikein_data[[col]] <= 0, na.rm = TRUE)) {
      warning(paste("列", colnames(log_observed_spikein_data)[col], "包含非正值，将在对数转换前处理"))
      # 将非正值替换为一个很小的正数，避免对数计算错误
      log_observed_spikein_data[[col]][log_observed_spikein_data[[col]] <= 0] <- 1e-10
    }
    
    log_observed_spikein_data[[col]] <- log2(log_observed_spikein_data[[col]])
  }
  
  # 保存log observed spike-in value.csv文件
  write.csv(log_observed_spikein_data, "log observed spike-in value.csv", row.names = FALSE)
  cat("对数转换后的观察值数据已保存为 'log observed spike-in value.csv'\n")
} else {
  # 如果没有数据，创建一个空的数据框
  write.csv(log_observed_spikein_data, "log observed spike-in value.csv", row.names = FALSE)
  cat("没有观察值数据可进行对数转换，创建空的 'log observed spike-in value.csv' 文件\n")
}

# 步骤4: 创建fpkm_matrix_endo.csv（剩余的数据）
# 如果没有找到任何匹配的基因，则整个fpkm_matrix都是剩余数据
if (length(extracted_genes) == 0) {
  fpkm_matrix_endo <- fpkm_matrix
  cat("没有找到任何匹配的基因，fpkm_matrix_endo.csv将包含所有原始数据\n")
} else {
  # 提取未被选中的行（剩余数据）
  fpkm_matrix_endo <- fpkm_matrix[-extracted_genes, ]
  cat(paste("提取了", length(extracted_genes), "个基因，剩余", nrow(fpkm_matrix_endo), "个基因\n"))
}

# 保存fpkm_matrix_endo.csv文件
write.csv(fpkm_matrix_endo, "fpkm_matrix_endo.csv", row.names = FALSE)

# 步骤5: 对fpkm_matrix_endo进行逐列筛选，保留大于0.1的数据
# 记录原始基因数量
original_count <- nrow(fpkm_matrix_endo)
cat(sprintf("\n开始筛选fpkm_matrix_endo数据，原始基因数量: %d\n", original_count))

# 确定样本列的范围（从第2列到最后一列）
sample_columns <- 2:ncol(fpkm_matrix_endo)

# 初始化筛选后的数据
fpkm_matrix_endo_filtered <- fpkm_matrix_endo

# 按顺序从第2列到最后一列依次筛选
for (i in seq_along(sample_columns)) {
  col_index <- sample_columns[i]
  col_name <- colnames(fpkm_matrix_endo_filtered)[col_index]
  
  # 筛选当前列大于0.1的行
  fpkm_matrix_endo_filtered <- fpkm_matrix_endo_filtered[fpkm_matrix_endo_filtered[, col_index] > 0.1, ]
  
  # 打印当前筛选信息
  cat(sprintf("筛选 %s 列后剩余基因数量: %d\n", 
              col_name, nrow(fpkm_matrix_endo_filtered)))
  
  # 如果没有基因剩余，停止筛选
  if (nrow(fpkm_matrix_endo_filtered) == 0) {
    cat("没有基因满足所有条件，停止筛选\n")
    break
  }
}

# 输出最终筛选后的数据
if (nrow(fpkm_matrix_endo_filtered) > 0) {
  write.csv(fpkm_matrix_endo_filtered, "fpkm_matrix_endo_filtered.csv", row.names = FALSE)
  
  # 计算筛选比例
  retention_rate <- round(nrow(fpkm_matrix_endo_filtered) / original_count * 100, 2)
  cat(sprintf("\n筛选完成，最终剩余基因数量: %d (保留率: %.2f%%)\n", 
              nrow(fpkm_matrix_endo_filtered), retention_rate))
  cat("筛选结果已保存为 'fpkm_matrix_endo_filtered.csv'\n")
} else {
  cat("没有基因满足所有筛选条件\n")
  # 如果没有基因满足条件，创建一个空的数据框
  fpkm_matrix_endo_filtered <- data.frame(matrix(ncol = ncol(fpkm_matrix_endo), nrow = 0))
  colnames(fpkm_matrix_endo_filtered) <- colnames(fpkm_matrix_endo)
  write.csv(fpkm_matrix_endo_filtered, "fpkm_matrix_endo_filtered.csv", row.names = FALSE)
}

# 步骤6: 创建log fpkm_matrix_endo_filtered.csv（对数转换的筛选后数据）
# 创建对数转换版本
log_fpkm_matrix_endo_filtered <- fpkm_matrix_endo_filtered

# 对数据列进行以2为底的对数转换（跳过第一列基因名）
if (nrow(log_fpkm_matrix_endo_filtered) > 0) {
  for (col in 2:ncol(log_fpkm_matrix_endo_filtered)) {
    # 检查是否有非正值，这些值不能取对数
    if (any(log_fpkm_matrix_endo_filtered[[col]] <= 0, na.rm = TRUE)) {
      warning(paste("列", colnames(log_fpkm_matrix_endo_filtered)[col], "包含非正值，将在对数转换前处理"))
      # 将非正值替换为一个很小的正数，避免对数计算错误
      log_fpkm_matrix_endo_filtered[[col]][log_fpkm_matrix_endo_filtered[[col]] <= 0] <- 1e-10
    }
    log_fpkm_matrix_endo_filtered[[col]] <- log2(log_fpkm_matrix_endo_filtered[[col]])
  }
  
  # 保存log fpkm_matrix_endo_filtered.csv文件
  write.csv(log_fpkm_matrix_endo_filtered, "log fpkm_matrix_endo_filtered.csv", row.names = FALSE)
  cat("对数转换后的数据已保存为 'log fpkm_matrix_endo_filtered.csv'\n")
} else {
  # 如果没有数据，创建一个空的数据框
  write.csv(log_fpkm_matrix_endo_filtered, "log fpkm_matrix_endo_filtered.csv", row.names = FALSE)
  cat("没有数据可进行对数转换，创建空的 'log fpkm_matrix_endo_filtered.csv' 文件\n")
}

# 输出信息
cat("\n第二部分处理完成！已生成以下文件:\n")
cat("1. input reference value.csv\n")
cat("2. log input reference value.csv\n")
cat("3. observed spike-in value.csv\n")
cat("4. log observed spike-in value.csv\n")
cat("5. fpkm_matrix_endo.csv\n")
cat("6. fpkm_matrix_endo_filtered.csv\n")
cat("7. log fpkm_matrix_endo_filtered.csv\n\n")
cat("fpkm_matrix文件总列数:", ncol(fpkm_matrix), "\n")
cat("fpkm_matrix文件数据列数:", data_cols_count, "\n")
cat("生成的文件列数:", ncol(new_ref_data), "\n")
cat("生成的文件行数:", nrow(new_ref_data), "\n")
cat("使用的reference列: B列和C列\n")
cat("乘法系数: 0.4\n")
cat("对数转换: 以2为底\n")
cat("从fpkm_matrix中提取的基因数:", length(extracted_genes), "\n")
cat("fpkm_matrix_endo中的基因数:", nrow(fpkm_matrix_endo), "\n")
cat("fpkm_matrix_endo_filtered中的基因数:", nrow(fpkm_matrix_endo_filtered), "\n")

# 显示一些统计信息
cat("\n原始数据统计 (前几列):\n")
print(summary(new_ref_data[, 2:min(5, ncol(new_ref_data))]))
cat("\n对数转换后数据统计 (前几列):\n")
print(summary(log_ref_data[, 2:min(5, ncol(log_ref_data))]))
cat("\n观察值数据统计 (前几列):\n")
print(summary(observed_spikein_data[, 2:min(5, ncol(observed_spikein_data))]))
cat("\n对数转换后观察值数据统计 (前几列):\n")
print(summary(log_observed_spikein_data[, 2:min(5, ncol(log_observed_spikein_data))]))
cat("\n剩余FPKM数据统计 (前几列):\n")
print(summary(fpkm_matrix_endo[, 2:min(5, ncol(fpkm_matrix_endo))]))
if (nrow(fpkm_matrix_endo_filtered) > 0) {
  cat("\n筛选后FPKM数据统计 (前几列):\n")
  print(summary(fpkm_matrix_endo_filtered[, 2:min(5, ncol(fpkm_matrix_endo_filtered))]))
  cat("\n对数转换后筛选FPKM数据统计 (前几列):\n")
  print(summary(log_fpkm_matrix_endo_filtered[, 2:min(5, ncol(log_fpkm_matrix_endo_filtered))]))
}

# ========== 第三部分：数据预处理、ME计算和线性拟合 ==========

cat("\n\n开始第三部分：数据预处理、ME计算和线性拟合\n")

# 读取数据文件（使用第二部分生成的文件）
observed_data <- read.csv("log observed spike-in value.csv", row.names = 1)
reference_data <- read.csv("log input reference value.csv", row.names = 1)

# 检查两个文件是否有相同的基因顺序
if(!identical(rownames(observed_data), rownames(reference_data))) {
  stop("两个文件的基因顺序不一致，请先确保基因顺序相同")
}

# 检查两个文件是否有相同的列数
if(ncol(observed_data) != ncol(reference_data)) {
  stop("两个文件的列数不一致")
}

# 计算每列数据的平均值
observed_means <- colMeans(observed_data, na.rm = TRUE)
reference_means <- colMeans(reference_data, na.rm = TRUE)

# 计算平均值的差值（reference - observed）
mean_differences <- reference_means - observed_means

# 将差值应用到observed数据中，实现均值校正
corrected_observed_data <- observed_data

# 对每一列加上对应的差值
for(i in 1:ncol(corrected_observed_data)) {
  corrected_observed_data[, i] <- corrected_observed_data[, i] + mean_differences[i]
}

# 验证校正后的均值是否与reference一致
corrected_means <- colMeans(corrected_observed_data, na.rm = TRUE)
verification <- round(corrected_means, 6) == round(reference_means, 6)

cat("校正验证结果（TRUE表示均值一致）：\n")
print(verification)

# 保存校正后的数据，使用指定的文件名
write.csv(corrected_observed_data, "log observed spike-in value_mean normalization.csv")

# 计算测量误差ME
# ME = (标准化后的观测值 - 参考值) / 参考值
ME_data <- (corrected_observed_data - reference_data) / reference_data

# 保存ME计算结果
write.csv(ME_data, "ME calculation.csv")

# 筛选ME数值在±0.05（含）之间的数据
# 获取基因名称列
gene_names <- rownames(ME_data)

# 创建一个工作簿来存储筛选和拟合结果
wb <- createWorkbook()

# 创建一个数据框来存储所有样本的拟合结果和基因数目
summary_results <- data.frame(
  Sample = character(),
  Gene_Count = integer(),
  Equation = character(),
  R_squared = numeric(),
  Remark = character(),  # 备注列
  stringsAsFactors = FALSE
)

# 对每个样本分别处理并添加到工作簿
for (col_name in colnames(ME_data)) {
  # 获取当前样本的ME数据
  current_me <- ME_data[[col_name]]
  
  # 筛选满足条件的基因和ME值（±0.05）
  valid_indices <- which(current_me >= -0.05 & current_me <= 0.05)
  valid_genes <- gene_names[valid_indices]
  valid_me <- current_me[valid_indices]
  
  # 从reference数据表中获取对应数据（x轴）
  reference_values <- numeric(length(valid_genes))
  for (i in seq_along(valid_genes)) {
    gene <- valid_genes[i]
    reference_values[i] <- reference_data[gene, col_name]
  }
  
  # 从原始observed数据表中获取对应数据（y轴）- 使用原始数据而非校正后的数据
  observed_values <- numeric(length(valid_genes))
  for (i in seq_along(valid_genes)) {
    gene <- valid_genes[i]
    observed_values[i] <- observed_data[gene, col_name]
  }
  
  # 创建当前样本的数据框
  sample_df <- data.frame(
    Gene = valid_genes,
    ME = valid_me,
    Reference_Value = reference_values,
    Observed_Value = observed_values  # 使用原始observed数据
  )
  
  # 进一步筛选：排除observed_value <= -6.64的数据点
  sample_df_filtered <- sample_df[sample_df$Observed_Value > -6.64, ]
  gene_count <- nrow(sample_df_filtered)
  
  # 添加工作表
  addWorksheet(wb, col_name)
  
  # 将数据写入工作表
  writeData(wb, sheet = col_name, x = sample_df_filtered)
  
  # 添加标题样式
  header_style <- createStyle(textDecoration = "bold", halign = "center")
  addStyle(wb, sheet = col_name, header_style, rows = 1, cols = 1:4, gridExpand = TRUE)
  
  # 添加边框
  border_style <- createStyle(border = "TopBottomLeftRight", borderColour = "#000000")
  addStyle(wb, sheet = col_name, border_style, rows = 1:(nrow(sample_df_filtered)+1), cols = 1:4, gridExpand = TRUE)
  
  # 进行线性拟合 (Observed_Value ~ Reference_Value)
  # 移除NA值
  valid_data <- sample_df_filtered[complete.cases(sample_df_filtered[, c("Reference_Value", "Observed_Value")]), ]
  valid_data_count <- nrow(valid_data)
  
  if (valid_data_count > 1) {
    # 执行线性回归
    fit <- lm(Observed_Value ~ Reference_Value, data = valid_data)
    
    # 提取拟合参数
    intercept <- coef(fit)[1]
    slope <- coef(fit)[2]
    r_squared <- summary(fit)$r.squared
    
    # 创建拟合方程字符串
    equation <- sprintf("y = %.4f + %.4f * x", intercept, slope)
    
    # 判断是否需要添加备注
    remark <- ""
    if (valid_data_count < 20 || r_squared < 0.95) {
      if (valid_data_count < 20 && r_squared < 0.95) {
        remark <- "数据点不足20且R²<0.95"
      } else if (valid_data_count < 20) {
        remark <- "数据点不足20"
      } else {
        remark <- "R²<0.95"
      }
    }
    
    # 将结果添加到汇总表
    summary_results <- rbind(summary_results, data.frame(
      Sample = col_name,
      Gene_Count = gene_count,
      Equation = equation,
      R_squared = r_squared,
      Remark = remark,
      stringsAsFactors = FALSE
    ))
    
    # 在工作表中添加拟合结果注释
    start_row <- nrow(sample_df_filtered) + 3
    writeData(wb, sheet = col_name, x = "样本信息:", startRow = start_row)
    writeData(wb, sheet = col_name, x = paste("样本名称:", col_name), startRow = start_row + 1)
    writeData(wb, sheet = col_name, x = paste("保留基因数目:", gene_count), startRow = start_row + 2)
    writeData(wb, sheet = col_name, x = paste("用于拟合的数据点数目:", valid_data_count), startRow = start_row + 3)
    writeData(wb, sheet = col_name, x = "线性拟合结果:", startRow = start_row + 5)
    writeData(wb, sheet = col_name, x = paste("模型: Observed_Value ~ Reference_Value"), startRow = start_row + 6)
    writeData(wb, sheet = col_name, x = paste("方程:", equation), startRow = start_row + 7)
    writeData(wb, sheet = col_name, x = paste("R²:", round(r_squared, 4)), startRow = start_row + 8)
    
    # 如果需要，在工作表中添加备注
    if (remark != "") {
      writeData(wb, sheet = col_name, x = paste("备注:", remark), startRow = start_row + 9)
    }
    
    cat(sprintf("样本 %s: 保留 %d 个基因, 拟合方程 %s, R² = %.4f\n", 
                col_name, gene_count, equation, r_squared))
  } else {
    cat(sprintf("样本 %s: 保留 %d 个基因, 但没有足够的数据进行线性拟合\n", col_name, gene_count))
    
    # 将结果添加到汇总表
    summary_results <- rbind(summary_results, data.frame(
      Sample = col_name,
      Gene_Count = gene_count,
      Equation = "无法拟合",
      R_squared = NA,
      Remark = "数据点不足，无法拟合",
      stringsAsFactors = FALSE
    ))
    
    # 在工作表中添加注释
    start_row <- nrow(sample_df_filtered) + 3
    writeData(wb, sheet = col_name, x = "样本信息:", startRow = start_row)
    writeData(wb, sheet = col_name, x = paste("样本名称:", col_name), startRow = start_row + 1)
    writeData(wb, sheet = col_name, x = paste("保留基因数目:", gene_count), startRow = start_row + 2)
    writeData(wb, sheet = col_name, x = "线性拟合结果:", startRow = start_row + 4)
    writeData(wb, sheet = col_name, x = "没有足够的数据进行线性拟合", startRow = start_row + 5)
  }
}

# 添加汇总工作表
addWorksheet(wb, "样本汇总")
writeData(wb, sheet = "样本汇总", x = summary_results)

# 设置汇总工作表的列宽和样式
setColWidths(wb, sheet = "样本汇总", cols = 1:5, widths = c(15, 15, 30, 15, 25))  # 增加第5列宽度
header_style <- createStyle(textDecoration = "bold", halign = "center")
addStyle(wb, sheet = "样本汇总", header_style, rows = 1, cols = 1:5, gridExpand = TRUE)
border_style <- createStyle(border = "TopBottomLeftRight", borderColour = "#000000")
addStyle(wb, sheet = "样本汇总", border_style, rows = 1:(nrow(summary_results)+1), cols = 1:5, gridExpand = TRUE)

# 保存工作簿，文件名为ME_filtered and linear fitting.xlsx
saveWorkbook(wb, "ME_filtered and linear fitting.xlsx", overwrite = TRUE)

# 同时保存筛选后的ME数据为CSV格式
# 创建一个数据框来存储所有样本筛选后的ME数据
ME_filtered_csv <- data.frame(Gene = gene_names)

for (col_name in colnames(ME_data)) {
  # 获取当前样本的ME数据
  current_me <- ME_data[[col_name]]
  
  # 筛选满足条件的基因和ME值（±0.05）
  valid_indices <- which(current_me >= -0.05 & current_me <= 0.05)
  valid_genes <- gene_names[valid_indices]
  valid_me <- current_me[valid_indices]
  
  # 创建一个临时向量，初始值为NA
  temp_col <- rep(NA, length(gene_names))
  names(temp_col) <- gene_names
  
  # 将符合条件的值填入
  temp_col[valid_genes] <- valid_me
  
  # 添加到数据框
  ME_filtered_csv[[col_name]] <- temp_col
}

# 保存筛选后的ME数据为CSV
write.csv(ME_filtered_csv, "ME calculation_filtered.csv", row.names = FALSE)

# 输出均值信息
cat("\n原始observed数据各列均值：\n")
print(observed_means)

cat("\nReference数据各列均值：\n")
print(reference_means)

cat("\n校正后的observed数据各列均值：\n")
print(corrected_means)

cat("\n各列应用的差值：\n")
print(mean_differences)

cat("\nME数据统计摘要：\n")
print(summary(as.matrix(ME_data)))

# 显示第三部分结果信息
cat("\n第三部分处理完成！\n")
cat("样本数量:", ncol(ME_data), "\n")
cat("校正后的数据已保存为: 'log observed spike-in value_mean normalization.csv'\n")
cat("测量误差ME数据已保存为: 'ME calculation.csv'\n")
cat("筛选后的ME数据已保存为: 'ME calculation_filtered.csv'\n")
cat("详细拟合结果已保存为: 'ME_filtered and linear fitting.xlsx'\n")

# 显示总体统计信息
cat("\n总体统计信息:\n")
cat(sprintf("平均每个样本保留基因数: %.1f\n", mean(summary_results$Gene_Count, na.rm = TRUE)))
if (nrow(summary_results) > 0) {
  cat(sprintf("最大保留基因数: %d (样本: %s)\n", 
              max(summary_results$Gene_Count, na.rm = TRUE), 
              summary_results$Sample[which.max(summary_results$Gene_Count)]))
  cat(sprintf("最小保留基因数: %d (样本: %s)\n", 
              min(summary_results$Gene_Count, na.rm = TRUE), 
              summary_results$Sample[which.min(summary_results$Gene_Count)]))
  
  # 显示需要备注的样本数量
  need_remark_count <- sum(summary_results$Remark != "", na.rm = TRUE)
  cat(sprintf("需要备注的样本数量: %d\n", need_remark_count))
}

# ========== 第四部分：计算spike-in基因拷贝数 ==========

cat("\n\n开始第四部分：计算spike-in基因拷贝数\n")

# 直接读取log observed spike-in value.csv文件（已经是log2转换后的数据）
log_spikein_data <- read.csv("log observed spike-in value.csv", header = TRUE)

# 获取基因名称列（第一列）
gene_names_spikein <- log_spikein_data[, 1]

# 获取所有样本的数据列（从第二列到最后一列）
sample_data_spikein <- log_spikein_data[, -1]

# 初始化拷贝数矩阵，与输入数据相同的结构
cp_matrix_spikein <- data.frame(Gene = gene_names_spikein)

# 对每个样本分别处理
for (i in 1:nrow(summary_results)) {
  sample_name <- summary_results$Sample[i]
  equation <- summary_results$Equation[i]
  
  # 检查是否有有效的拟合方程且样本存在于spike-in基因数据中
  if (equation != "无法拟合" && sample_name %in% colnames(sample_data_spikein)) {
    # 从方程字符串中提取斜率和截距
    # 方程格式应为: "y = [截距] + [斜率] * x"
    parts <- strsplit(equation, " ")[[1]]
    intercept <- as.numeric(parts[3])
    slope <- as.numeric(parts[5])
    
    # 获取当前样本的y值（已经是log2(FPKM)值）
    y_values <- sample_data_spikein[[sample_name]]
    
    # 初始化拷贝数值向量
    cp_values <- numeric(length(y_values))
    
    # 对每个基因的值分别处理
    for (j in 1:length(y_values)) {
      y_val <- y_values[j]
      
      # 检查是否小于等于-6.64
      if (!is.na(y_val) && y_val <= -6.64) {
        # 如果小于等于-6.64，拷贝数设为0
        cp_values[j] <- 0
      } else if (!is.na(y_val)) {
        # 否则使用拟合方程计算x值 (x = (y - 截距) / 斜率)
        # x是对数拷贝数
        logCP_value <- (y_val - intercept) / slope
        
        # 计算拷贝数 (2^x)
        cp_values[j] <- 2^logCP_value
      } else {
        # 如果是NA值，保持为NA
        cp_values[j] <- NA
      }
    }
    
    # 将拷贝数格式化为科学计数法，保留2位小数
    # 对于0值，直接显示为0.00e+00
    cp_values_formatted <- ifelse(cp_values == 0, "0.00e+00", formatC(cp_values, format = "e", digits = 2))
    
    # 将结果添加到拷贝数矩阵
    cp_matrix_spikein[[sample_name]] <- cp_values_formatted
    
    cat(sprintf("样本 %s: 已计算spike-in基因拷贝数值\n", sample_name))
  } else {
    # 如果没有有效的拟合方程，填充NA值
    cp_matrix_spikein[[sample_name]] <- rep(NA, length(gene_names_spikein))
    if (equation == "无法拟合") {
      cat(sprintf("样本 %s: 没有有效的拟合方程，填充NA值\n", sample_name))
    } else if (!(sample_name %in% colnames(sample_data_spikein))) {
      cat(sprintf("样本 %s: 在spike-in基因数据中不存在，填充NA值\n", sample_name))
    }
  }
}

# 保存spike-in拷贝数矩阵为CSV文件
write.csv(cp_matrix_spikein, "cp_matrix_spike-in.csv", row.names = FALSE)

# 显示第四部分结果信息
cat("\n第四部分处理完成！\n")
cat("spike-in基因数量:", nrow(cp_matrix_spikein), "\n")
cat("处理的样本数量:", sum(summary_results$Equation != "无法拟合" & summary_results$Sample %in% colnames(sample_data_spikein)), "\n")
cat("跳过的样本数量:", sum(summary_results$Equation == "无法拟合" | !(summary_results$Sample %in% colnames(sample_data_spikein))), "\n")
cat("文件已保存为 cp_matrix_spike-in.csv\n")

# 显示矩阵的前几行和前几列
cat("\nspike-in拷贝数矩阵预览 (前5行，前5列):\n")
print(head(cp_matrix_spikein[, 1:min(6, ncol(cp_matrix_spikein))], 5))

# ========== 第五部分：计算内源基因拷贝数 ==========

cat("\n\n开始第五部分：计算内源基因拷贝数\n")

# 读取内源基因FPKM数据（使用第二部分生成的文件）
endo_fpkm_data <- read.csv("log fpkm_matrix_endo_filtered.csv", header = TRUE)

# 获取基因名称列（第一列）
gene_names_endo <- endo_fpkm_data[, 1]

# 获取所有样本的数据列（从第二列到最后一列）
sample_data_endo <- endo_fpkm_data[, -1]

# 初始化拷贝数矩阵，与输入数据相同的结构
cp_matrix <- data.frame(Gene = gene_names_endo)

# 对每个样本分别处理
for (i in 1:nrow(summary_results)) {
  sample_name <- summary_results$Sample[i]
  equation <- summary_results$Equation[i]
  
  # 检查是否有有效的拟合方程且样本存在于内源基因数据中
  if (equation != "无法拟合" && sample_name %in% colnames(sample_data_endo)) {
    # 从方程字符串中提取斜率和截距
    # 方程格式应为: "y = [截距] + [斜率] * x"
    parts <- strsplit(equation, " ")[[1]]
    intercept <- as.numeric(parts[3])
    slope <- as.numeric(parts[5])
    
    # 获取当前样本的y值（内源基因的log FPKM值）
    y_values <- sample_data_endo[[sample_name]]
    
    # 使用拟合方程计算x值 (x = (y - 截距) / 斜率)
    # x是对数拷贝数
    logCP_values <- (y_values - intercept) / slope
    
    # 计算拷贝数 (2^x)
    cp_values <- 2^logCP_values
    
    # 将拷贝数格式化为科学计数法，保留2位小数
    cp_values_formatted <- formatC(cp_values, format = "e", digits = 2)
    
    # 将结果添加到拷贝数矩阵
    cp_matrix[[sample_name]] <- cp_values_formatted
    
    cat(sprintf("样本 %s: 已计算内源基因拷贝数值\n", sample_name))
  } else {
    # 如果没有有效的拟合方程，填充NA值
    cp_matrix[[sample_name]] <- rep(NA, length(gene_names_endo))
    if (equation == "无法拟合") {
      cat(sprintf("样本 %s: 没有有效的拟合方程，填充NA值\n", sample_name))
    } else if (!(sample_name %in% colnames(sample_data_endo))) {
      cat(sprintf("样本 %s: 在内源基因数据中不存在，填充NA值\n", sample_name))
    }
  }
}

# 保存内源基因拷贝数矩阵为CSV文件
write.csv(cp_matrix, "cp_matrix_endo.csv", row.names = FALSE)

# 显示第五部分结果信息
cat("\n第五部分处理完成！\n")
cat("内源基因数量:", nrow(cp_matrix), "\n")
cat("处理的样本数量:", sum(summary_results$Equation != "无法拟合" & summary_results$Sample %in% colnames(sample_data_endo)), "\n")
cat("跳过的样本数量:", sum(summary_results$Equation == "无法拟合" | !(summary_results$Sample %in% colnames(sample_data_endo))), "\n")
cat("文件已保存为 cp_matrix_endo.csv\n")

# 显示矩阵的前几行和前几列
cat("\n内源基因拷贝数矩阵预览 (前5行，前5列):\n")
print(head(cp_matrix[, 1:min(6, ncol(cp_matrix))], 5))

# ========== 第六部分：计算拷贝数比值 ==========

cat("\n\n开始第六部分：计算拷贝数比值\n")

# 功能1：计算拷贝数比值（cp_matrix_spike-in.csv与input reference value.csv的比值）

# 读取spike-in拷贝数矩阵
cp_spikein_data <- read.csv("cp_matrix_spike-in.csv", header = TRUE)

# 读取input reference value.csv文件
input_ref_data <- read.csv("input reference value.csv", header = TRUE)

# 获取基因名称列（第一列）
gene_names_cp <- cp_spikein_data[, 1]

# 初始化拷贝数比值矩阵
cp_ratio_matrix <- data.frame(Gene = gene_names_cp)

# 对每个样本分别处理
for (col_name in colnames(cp_spikein_data)[-1]) {
  if (col_name %in% colnames(input_ref_data)) {
    # 获取当前样本的spike-in拷贝数
    cp_values <- cp_spikein_data[[col_name]]
    
    # 将科学计数法字符串转换为数值
    cp_values_numeric <- as.numeric(cp_values)
    
    # 获取当前样本的输入参考值
    ref_values <- input_ref_data[[col_name]]
    
    # 计算拷贝数比值（spike-in拷贝数 / 输入参考值）
    # 注意：如果spike-in拷贝数为0，比值为0；如果输入参考值为0，比值为NA
    ratio_values <- numeric(length(cp_values_numeric))
    
    for (i in 1:length(cp_values_numeric)) {
      cp_val <- cp_values_numeric[i]
      ref_val <- ref_values[i]
      
      if (!is.na(cp_val) && !is.na(ref_val)) {
        if (cp_val == 0) {
          ratio_values[i] <- 0
        } else if (ref_val == 0) {
          ratio_values[i] <- NA
        } else {
          ratio_values[i] <- cp_val / ref_val
        }
      } else {
        ratio_values[i] <- NA
      }
    }
    
    # 将比值保留2位小数
    ratio_values_formatted <- ifelse(ratio_values == 0, "0.00", 
                                     formatC(ratio_values, format = "f", digits = 2))
    
    # 将结果添加到比值矩阵
    cp_ratio_matrix[[col_name]] <- ratio_values_formatted
    
    cat(sprintf("样本 %s: 已计算拷贝数比值\n", col_name))
  } else {
    # 如果样本在input reference中不存在，填充NA值
    cp_ratio_matrix[[col_name]] <- rep(NA, length(gene_names_cp))
    cat(sprintf("样本 %s: 在input reference中不存在，填充NA值\n", col_name))
  }
}

# 保存拷贝数比值矩阵为CSV文件
write.csv(cp_ratio_matrix, "cp ratio_matrix_calibrated and input spike-in.csv", row.names = FALSE)

cat("\n拷贝数比值矩阵已保存为 'cp ratio_matrix_calibrated and input spike-in.csv'\n")

# 功能2：在log水平计算比值

# 步骤1：生成log cp_matrix_spike-in.csv
cat("\n生成log cp_matrix_spike-in.csv...\n")

# 初始化log拷贝数矩阵
log_cp_matrix_spikein <- data.frame(Gene = gene_names_cp)

# 对每个样本分别处理
for (col_name in colnames(cp_spikein_data)[-1]) {
  # 获取当前样本的spike-in拷贝数
  cp_values <- cp_spikein_data[[col_name]]
  
  # 将科学计数法字符串转换为数值
  cp_values_numeric <- as.numeric(cp_values)
  
  # 计算log2拷贝数
  # 如果拷贝数为0，log值仍记录为0
  log_cp_values <- numeric(length(cp_values_numeric))
  
  for (i in 1:length(cp_values_numeric)) {
    cp_val <- cp_values_numeric[i]
    
    if (!is.na(cp_val)) {
      if (cp_val == 0) {
        log_cp_values[i] <- 0  # 拷贝数为0时，log值记录为0
      } else {
        log_cp_values[i] <- log2(cp_val)
      }
    } else {
      log_cp_values[i] <- NA
    }
  }
  
  # 将log拷贝数添加到矩阵
  log_cp_matrix_spikein[[col_name]] <- log_cp_values
  
  cat(sprintf("样本 %s: 已计算log拷贝数值\n", col_name))
}

# 保存log拷贝数矩阵为CSV文件
write.csv(log_cp_matrix_spikein, "log cp_matrix_spike-in.csv", row.names = FALSE)

cat("log拷贝数矩阵已保存为 'log cp_matrix_spike-in.csv'\n")

# 步骤2：计算log水平的拷贝数比值
cat("\n计算log水平的拷贝数比值...\n")

# 读取log input reference value.csv文件
log_input_ref_data <- read.csv("log input reference value.csv", header = TRUE)

# 初始化log拷贝数比值矩阵
log_cp_ratio_matrix <- data.frame(Gene = gene_names_cp)

# 对每个样本分别处理
for (col_name in colnames(log_cp_matrix_spikein)[-1]) {
  if (col_name %in% colnames(log_input_ref_data)) {
    # 获取当前样本的log spike-in拷贝数
    log_cp_values <- log_cp_matrix_spikein[[col_name]]
    
    # 获取当前样本的log输入参考值
    log_ref_values <- log_input_ref_data[[col_name]]
    
    # 计算log水平的拷贝数比值 (log2(校准后拷贝数)) / (log2(输入拷贝数))
    log_ratio_values <- numeric(length(log_cp_values))
    
    for (i in 1:length(log_cp_values)) {
      log_cp_val <- log_cp_values[i]
      log_ref_val <- log_ref_values[i]
      
      if (!is.na(log_cp_val) && !is.na(log_ref_val)) {
        # 如果log spike-in拷贝数为0（原始拷贝数为0）或log输入参考值为0，比值为NA
        if (log_cp_val == 0 || log_ref_val == 0) {
          log_ratio_values[i] <- NA
        } else {
          log_ratio_values[i] <- log_cp_val / log_ref_val
        }
      } else {
        log_ratio_values[i] <- NA
      }
    }
    
    # 将log比值添加到矩阵
    log_cp_ratio_matrix[[col_name]] <- log_ratio_values
    
    cat(sprintf("样本 %s: 已计算log拷贝数比值\n", col_name))
  } else {
    # 如果样本在log input reference中不存在，填充NA值
    log_cp_ratio_matrix[[col_name]] <- rep(NA, length(gene_names_cp))
    cat(sprintf("样本 %s: 在log input reference中不存在，填充NA值\n", col_name))
  }
}

# 保存log拷贝数比值矩阵为CSV文件
write.csv(log_cp_ratio_matrix, "log cp ratio_matrix_calibrated and input spike-in.csv", row.names = FALSE)

cat("log拷贝数比值矩阵已保存为 'log cp ratio_matrix_calibrated and input spike-in.csv'\n")

# 显示第六部分结果信息
cat("\n第六部分处理完成！\n")
cat("基因数量:", nrow(cp_ratio_matrix), "\n")
cat("处理的样本数量:", sum(colnames(cp_spikein_data)[-1] %in% colnames(input_ref_data)), "\n")
cat("生成的比值矩阵文件:\n")
cat("1. cp ratio_matrix_calibrated and input spike-in.csv\n")
cat("2. log cp_matrix_spike-in.csv\n")
cat("3. log cp ratio_matrix_calibrated and input spike-in.csv\n")

# 显示矩阵的前几行和前几列
cat("\n拷贝数比值矩阵预览 (前5行，前5列):\n")
print(head(cp_ratio_matrix[, 1:min(6, ncol(cp_ratio_matrix))], 5))

cat("\nlog拷贝数矩阵预览 (前5行，前5列):\n")
print(head(log_cp_matrix_spikein[, 1:min(6, ncol(log_cp_matrix_spikein))], 5))

cat("\nlog拷贝数比值矩阵预览 (前5行，前5列):\n")
print(head(log_cp_ratio_matrix[, 1:min(6, ncol(log_cp_ratio_matrix))], 5))

# 总体完成信息
cat("\n所有处理已完成！\n")
cat("整个流程从count数据到拷贝数计算和比值分析已全部完成！\n")
cat("==========================================\n")
cat("总结：生成了以下重要文件：\n")
cat("1. gene_lengths.csv - 基因长度信息\n")
cat("2. fpkm_matrix.csv - FPKM表达矩阵\n")
cat("3. log2_fpkm_matrix.csv - log2转换的FPKM矩阵\n")
cat("4. input reference value.csv - 输入参考值\n")
cat("5. log input reference value.csv - 对数转换的输入参考值\n")
cat("6. observed spike-in value.csv - 观察到的spike-in值\n")
cat("7. log observed spike-in value.csv - 对数转换的观察值\n")
cat("8. fpkm_matrix_endo.csv - 内源基因FPKM矩阵\n")
cat("9. fpkm_matrix_endo_filtered.csv - 筛选后的内源基因FPKM矩阵\n")
cat("10. log fpkm_matrix_endo_filtered.csv - 对数转换的筛选后FPKM矩阵\n")
cat("11. log observed spike-in value_mean normalization.csv - 均值标准化后的观察值\n")
cat("12. ME calculation.csv - 测量误差数据\n")
cat("13. ME calculation_filtered.csv - 筛选后的ME数据\n")
cat("14. ME_filtered and linear fitting.xlsx - 详细的拟合结果\n")
cat("15. cp_matrix_spike-in.csv - spike-in基因拷贝数矩阵\n")
cat("16. cp_matrix_endo.csv - 内源基因拷贝数矩阵\n")
cat("17. cp ratio_matrix_calibrated and input spike-in.csv - 拷贝数比值矩阵\n")
cat("18. log cp_matrix_spike-in.csv - log拷贝数矩阵\n")
cat("19. log cp ratio_matrix_calibrated and input spike-in.csv - log拷贝数比值矩阵\n")
cat("==========================================\n")