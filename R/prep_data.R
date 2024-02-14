prep_rna_data = function(rna_data,
												 sample_info)
{
	# tar_load(rna_data)
	# rna_data = rna_data[, seq_len(16)]
	# tar_load(sample_info)
	rna_matrix = as.matrix(rna_data[, -1])
	rownames(rna_matrix) = rna_data$gene_id
	sample_info = sample_info |>
		dplyr::filter(sample_id %in% colnames(rna_matrix))
	rna_matrix = rna_matrix[, sample_info$sample_id]
	
	sample_info$treatment = factor(sample_info$treatment, levels = c("normal_adjacent", "cancerous"))
	rna_dds = DESeq2::DESeqDataSetFromMatrix(rna_matrix, sample_info,
																					 design = ~ patient + treatment)
	rna_dds = DESeq2::estimateSizeFactors(rna_dds)
	rna_dds
}