rename_experimental_samples = function(metabolomics_dataset)
{
	# metabolomics_dataset = tar_read(bioamines)
	names(metabolomics_dataset) = gsub("^x", "s", names(metabolomics_dataset))
	metabolomics_dataset
}

median_normalization = function(feature_abundance)
{
	# feature_abundance = tar_read(lipidomics)$abundance
	long_abundance = feature_abundance |>
		tidyr::pivot_longer(cols = -feature_id,
												names_to = "sample_id",
												values_to = "abundance")
	median_values = long_abundance |>
		dplyr::group_by(sample_id) |>
		dplyr::summarise(median = median(abundance, na.rm = TRUE))
	
	tmp_abundance = dplyr::left_join(long_abundance, median_values, by = c("sample_id"))
	out_abundance = tmp_abundance |>
		dplyr::mutate(abundance = abundance / median,
									abundance_type = "median_normalized") |>
		dplyr::select(-median)
	out_abundance
}

matrix_2_long = function(wide_matrix)
{
	# wide_matrix = tar_load(rna_norm)
	tmp_df = tibble::as_tibble(wide_matrix) |>
		dplyr::mutate(feature_id = rownames(wide_matrix))
	long_df = tmp_df |>
		tidyr::pivot_longer(cols = tidyselect::where(is.numeric), names_to = "sample_id",
												values_to = "abundance")
	long_df
}

long_2_matrix = function(long_data)
{
	# long_data = tar_read(rna_norm)
	wide_df = long_data |>
		tidyr::pivot_wider(id_cols = feature_id,
											 names_from = "sample_id",
											 values_from = "abundance")
	wide_matrix = as.matrix(wide_df |> dplyr::select(-feature_id))
	rownames(wide_matrix) = wide_df$feature_id
	wide_matrix
}

keep_presence = function(normalized_data,
												 sample_info,
												 fraction = 0.25)
{
	# normalized_data = tar_read(bioamines_norm)
	# tar_load(sample_info)
	# fraction = 0.25
	norm_wide = long_2_matrix(normalized_data)
	sample_info = sample_info |>
		dplyr::filter(sample_id %in% colnames(norm_wide))
	norm_wide = norm_wide[, sample_info$sample_id]
	
	keep_norm = visualizationQualityControl::keep_non_missing_percentage(norm_wide, sample_info$treatment,
																																			 keep_num = fraction,
																																			 missing_value = c(NA, 0))
	keep_features = names(keep_norm)[keep_norm]
	normalized_out = normalized_data |>
		dplyr::filter(feature_id %in% keep_features)
	return(normalized_out)
}

sample_correlations = function(normalized_data)
{
	# normalized_data = tar_read(bioamines_keep)
	wide_matrix = long_2_matrix(normalized_data)
	
	sample_sample_cor = ICIKendallTau::ici_kendalltau(wide_matrix, global_na = c(0, NA), perspective = "global", scale_max = TRUE, diag_good = TRUE)
	
	if (any(is.na(wide_matrix))) {
		wide_matrix[is.na(wide_matrix)] = 0
	}
	
	sample_all_pca = prcomp(t(log1p(wide_matrix)), center = TRUE, scale. = FALSE)
	pool_blank_samples = grepl("^pool|^blank", colnames(wide_matrix))
	sample_noblanks_pca = prcomp(t(log1p(wide_matrix[, !pool_blank_samples])), center = TRUE, scale. = FALSE)
	
	return(list(correlation = sample_sample_cor,
							all_pca = sample_all_pca,
							noblanks_pca = sample_noblanks_pca))
}