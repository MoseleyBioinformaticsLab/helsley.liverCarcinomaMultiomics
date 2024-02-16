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

split_intensities_feature_metadata = function(all_data)
{
	# tar_load(lipidomics_file)
	# all_data = suppressMessages(readxl::read_excel(lipidomics_file,
	# sheet = "Data",
	# skip = 8)) |>
	# janitor::clean_names()  |>
	# rename_experimental_samples()
	start_samples = which(grepl("^s[[:digit:]]", colnames(all_data)))[1]
	sample_locs = seq(start_samples, ncol(all_data))
	metadata_locs = seq_len(start_samples - 1)
	
	if ("bin_base_name" %in% colnames(all_data)) {
		feature_id = janitor::make_clean_names(all_data[["bin_base_name"]])
	} else {
		feature_id = janitor::make_clean_names(all_data[["identifier"]])
	}
	
	sample_data = all_data[, sample_locs]
	sample_data = purrr::map(sample_data, as.numeric) |>
		dplyr::bind_cols()
	metadata = all_data[, metadata_locs]
	sample_data$feature_id = feature_id
	metadata$feature_id = feature_id
	
	list(abundance = sample_data,
			 metadata = metadata)
}