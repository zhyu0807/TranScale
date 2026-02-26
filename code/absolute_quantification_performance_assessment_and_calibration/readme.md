This code operates within the R environment to perform absolute quantification performance assessment and data calibration for transcriptomics sequencing, based on the standard values of TranScale reference materials. 
The workflow ultimately generates an absolute copy number expression matrix.

Input Files:
TranScale Standard Reference Values: A data table containing the certified reference values for the calibrators.
Transcriptome Count Matrix: Raw count data generated from transcriptomics sequencing.

Output Files:
gene_lengths.csv: Gene length information.
fpkm_matrix.csv: FPKM (Fragments Per Kilobase of transcript per Million mapped reads) expression matrix.
log2_fpkm_matrix.csv: $log_2$-transformed FPKM matrix.
input reference value.csv: Input reference values for the calibrators.
log input reference value.csv: Log-transformed input reference values.
observed spike-in value.csv: Observed values for the spike-in calibrators.
log observed spike-in value.csv: Log-transformed observed spike-in values.
fpkm_matrix_endo.csv: FPKM matrix for endogenous genes.
fpkm_matrix_endo_filtered.csv: Filtered FPKM matrix for endogenous genes.
log fpkm_matrix_endo_filtered.csv: Log-transformed filtered FPKM matrix.
log observed spike-in value_mean normalization.csv: Mean-normalized observed values for spike-ins.
ME calculation.csv: Measurement Error (ME) calculation results.
ME calculation_filtered.csv: Filtered Measurement Error (ME) data.
ME_filtered and linear fitting.xlsx: Detailed linear regression and fitting results.
cp_matrix_spike-in.csv: Absolute copy number matrix for spike-in genes.
cp_matrix_endo.csv: Absolute copy number matrix for endogenous genes.
cp ratio_matrix_calibrated and input spike-in.csv: Matrix of calibrated copy number ratios.
log cp_matrix_spike-in.csv: Log-transformed copy number matrix for spike-ins.
log cp ratio_matrix_calibrated and input spike-in.csv: Log-transformed copy number ratio matrix.
