prep_rna_data_patient = function(rna_data,
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

prep_rna_data_treatment = function(rna_data,
																	 sample_info)
{
	rna_matrix = as.matrix(rna_data[, -1])
	rownames(rna_matrix) = rna_data$gene_id
	sample_info = sample_info |>
		dplyr::filter(sample_id %in% colnames(rna_matrix))
	rna_matrix = rna_matrix[, sample_info$sample_id]
	
	sample_info$treatment = factor(sample_info$treatment, levels = c("normal_adjacent", "cancerous"))
	rna_dds = DESeq2::DESeqDataSetFromMatrix(rna_matrix, sample_info,
																					 design = ~ treatment)
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

determine_matched_samples = function(sample_info)
{
	# tar_load(sample_info)
	# sample_info = sample_info |> dplyr::filter(!sample_id %in% "s11")
	n_sample_patient = sample_info |>
		dplyr::group_by(patient) |>
		dplyr::summarise(n_sample = dplyr::n()) |>
		dplyr::filter(n_sample == 2)
	sample_info_out = sample_info |>
		dplyr::filter(patient %in% n_sample_patient$patient)
	sample_info_out
}

calculate_patient_ratios = function(normalized_data,
																		sample_info)
{
	# normalized_data = tar_read(rna_keep)
	# tar_load(sample_info)
	# sample_info = sample_info |> dplyr::filter(!sample_id %in% "s11")
	intersect_samples = base::intersect(unique(normalized_data$sample_id), sample_info$sample_id)
	matched_sample_info = determine_matched_samples(sample_info |> 
																										dplyr::filter(sample_id %in% intersect_samples))
	normalized_data = dplyr::inner_join(normalized_data, matched_sample_info[, c("patient", "sample_id", "treatment")],
																			by = "sample_id")
	patient_ratios = normalized_data |>
		tidyr::pivot_wider(id_cols = c("patient", "feature_id"),
											 names_from = "treatment",
											 values_from = "abundance") |>
		dplyr::mutate(abundance = cancerous / normal_adjacent,
									sample_id = patient)
	patient_ratios
}