prep_rna_data_patient = function(rna_data = rna_data, sample_metadata = sample_info,
																 counts_locs = seq(2, 16),
																 feature_metadata = seq(17, 25))
{
	# tar_load(rna_data)
	# sample_metadata = tar_read(sample_info)
	# count_locs = seq(2, 16)
	# feature_metadata = seq(17, 25)
	rna_matrix = as.matrix(rna_data[, count_locs])
	rownames(rna_matrix) = rna_data$gene_id
	metadata_cols = c("gene_id", names(rna_data)[feature_metadata])
	gene_metadata = rna_data[, metadata_cols]
	sample_metadata = sample_metadata |>
		dplyr::filter(sample_id %in% colnames(rna_matrix))
	sample_metadata_matched = determine_matched_samples(sample_metadata)
	rna_matrix = rna_matrix[, sample_metadata_matched$sample_id]
	
	gene_ranges = GRanges(seqnames = gene_metadata$gene_chr,
												ranges = IRanges(start = gene_metadata$gene_start,
																				 end = gene_metadata$gene_end),
												strand = gene_metadata$gene_strand,
												feature_id = gene_metadata$gene_id,
												biotype = gene_metadata$gene_biotype,
												description = gene_metadata$gene_description)
	
	sample_metadata$treatment = factor(sample_metadata$treatment, levels = c("normal_adjacent", "cancerous"))
	sample_metadata$patient = factor(sample_metadata$patient)
	rna_dds = DESeq2::DESeqDataSetFromMatrix(rna_matrix, sample_metadata,
																					 design = ~ patient + treatment,
																					 rowRanges = gene_ranges)
	rna_dds = DESeq2::estimateSizeFactors(rna_dds)
	rna_dds
}

prep_rna_data_treatment = function(rna_data = rna_data, sample_metadata = sample_info,
																	 count_locs = seq(2, 16),
																	 feature_metadata = seq(17, 25))
{
	rna_matrix = as.matrix(rna_data[, count_locs])
	rownames(rna_matrix) = rna_data$gene_id
	metadata_cols = c("gene_id", names(rna_data)[feature_metadata])
	gene_metadata = rna_data[, metadata_cols]
	sample_metadata = sample_metadata |>
		dplyr::filter(sample_id %in% colnames(rna_matrix))
	rna_matrix = rna_matrix[, sample_metadata$sample_id]
	
	# anything below min_count is set to 0
	
	gene_ranges = GRanges(seqnames = gene_metadata$gene_chr,
												ranges = IRanges(start = gene_metadata$gene_start,
																				 end = gene_metadata$gene_end),
												strand = gene_metadata$gene_strand,
												feature_id = gene_metadata$gene_id,
												name = gene_metadata$gene_name,
												biotype = gene_metadata$gene_biotype,
												description = gene_metadata$gene_description)
	
	sample_metadata$treatment = factor(sample_metadata$treatment, levels = c("normal_adjacent", "cancerous"))
	sample_metadata$patient = factor(sample_metadata$patient)
	rna_dds = DESeq2::DESeqDataSetFromMatrix(rna_matrix, sample_metadata,
																					 design = ~ treatment,
																					 rowRanges = gene_ranges)
	rna_dds = DESeq2::estimateSizeFactors(rna_dds)
	rna_dds
}

setup_metabolomics = function(all_data, 
															sample_info,
															metabolite_type = "",
															replace_zero = TRUE)
{
	# tar_load(lipidomics_file)
	# all_data = suppressMessages(readxl::read_excel(lipidomics_file,
	# sheet = "Data",
	# skip = 8)) |>
	# janitor::clean_names()  |>
	# rename_experimental_samples()
	# tar_load(sample_info)
	# metabolite_type = "lipidomics"
	start_samples = which(grepl("^s[[:digit:]]", colnames(all_data)))[1]
	sample_locs = seq(start_samples, ncol(all_data))
	metadata_locs = seq_len(start_samples - 1)
	
	if ("bin_base_name" %in% colnames(all_data)) {
		feature_id = paste0(janitor::make_clean_names(all_data[["bin_base_name"]]), ".",
												metabolite_type)
	} else {
		feature_id = paste0(janitor::make_clean_names(all_data[["identifier"]]), ".",
												metabolite_type)
	}
	
	if ("in_ch_i_key" %in% colnames(all_data)) {
		all_data = dplyr::rename(all_data, "in_chi_key" = "in_ch_i_key")
		message("renamed inchi keys")
	}
	sample_data = all_data[, sample_locs]
	sample_data = purrr::map(sample_data, as.numeric) |>
		dplyr::bind_cols()
	metadata = all_data[, metadata_locs]
	sample_data$feature_id = feature_id
	
	sample_matrix = sample_data |>
		dplyr::select(-feature_id) |>
		as.matrix()
	rownames(sample_matrix) = sample_data$feature_id
	
	sample_info_extra = add_blanks_pooled(colnames(sample_matrix), sample_info)
	sample_matrix = sample_matrix[, sample_info_extra$sample_id]
	
	sample_matrix[is.na(sample_matrix)] = 0
	
	metadata$feature_id = feature_id
	rownames(metadata) = metadata$feature_id
	
	out_data = DESeq2::DESeqDataSetFromMatrix(sample_matrix, sample_info_extra, design = ~ treatment)
	
	metadata_df = DataFrame(metadata)
	mcols(out_data) = metadata_df
	
	out_data = DESeq2::estimateSizeFactors(out_data)
	out_data
}

add_blanks_pooled = function(sample_names,
														 sample_info)
{
	# sample_names = colnames(cor_pca$correlation$cor)
	extra_samples = tibble::tibble(sample_id = setdiff(sample_names, sample_info$sample_id)) |>
		dplyr::mutate(treatment = dplyr::case_when(
			grepl("^pool", sample_id) ~ "pooled",
			grepl("^blank", sample_id) ~ "blank"
		),
		patient = "none")
	sample_info = dplyr::bind_rows(sample_info,
																 extra_samples)
	sample_info
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
