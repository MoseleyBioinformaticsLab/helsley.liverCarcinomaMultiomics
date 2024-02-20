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
												 fraction = 0.75)
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

floor_values = function(rna_matrix)
{
	new_matrix = rna_matrix
	new_matrix[rna_matrix < 10] = 0
	new_matrix
}

keep_presence_dds = function(rna_dds,
												 fraction = 0.75)
{
	# rna_dds = tar_read(rna_dds_treatment)
	# fraction = 0.75
	norm_wide = counts(rna_dds)
	sample_info = as.data.frame(colData(rna_dds))
	# set anything less than 10 to 0 so it still counts as missing
	norm_wide[norm_wide < 10] = 0
	keep_norm = visualizationQualityControl::keep_non_missing_percentage(norm_wide, sample_info$treatment,
																																			 keep_num = fraction,
																																			 missing_value = c(NA, 0))
	rna_dds = rna_dds[keep_norm, ]
	return(rna_dds)
}

sample_correlations_pca = function(normalized_data,
																	 sample_info)
{
	# normalized_data = tar_read(bioamines_keep)
	# tar_load(sample_info)
	# 
	# normalized_data = tar_read(rna_ratios)
	# sample_info = tar_read(patient_info)
	wide_matrix = long_2_matrix(normalized_data)
	wide_matrix[is.nan(wide_matrix) | is.infinite(wide_matrix)] = NA
	
	sample_sample_cor = ICIKendallTau::ici_kendalltau(wide_matrix, global_na = c(0, NA), perspective = "global", scale_max = TRUE, diag_good = TRUE)
	
	if (any(is.na(wide_matrix))) {
		wide_matrix[is.na(wide_matrix)] = 0
	}
	
	sample_all_pca = prcomp(t(log1p(wide_matrix)), center = TRUE, scale. = FALSE)
	pool_blank_samples = grepl("^pool|^blank", colnames(wide_matrix))
	sample_noblanks_pca = prcomp(t(log1p(wide_matrix[, !pool_blank_samples])), center = TRUE, scale. = FALSE)
	
	match_samples = intersect(colnames(wide_matrix), sample_info$sample_id)
	sample_info2 = sample_info |>
		dplyr::filter(sample_id %in% match_samples)
	matrix_cor = sample_sample_cor$cor[sample_info2$sample_id, sample_info2$sample_id]
	median_cor = visualizationQualityControl::median_correlations(matrix_cor, sample_info2$treatment)
	median_outliers = visualizationQualityControl::determine_outliers(median_correlations = median_cor)
	use_samples = median_outliers |>
		dplyr::filter(!outlier) |>
		dplyr::pull(sample_id)
	
	sample_nooutlier_pca = prcomp(t(log1p(wide_matrix[, use_samples])), center = TRUE, scale. = FALSE)
	
	return(list(correlation = sample_sample_cor,
							all_pca = sample_all_pca,
							noblanks_pca = sample_noblanks_pca,
							nooutlier_pca = sample_nooutlier_pca))
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

create_qcqa_plots = function(cor_pca,
														 sample_info,
														 color_scales)
{
	# cor_pca = tar_read(bioamines_cor_pca)
	# cor_pca = tar_read(primary_metabolism_cor_pca)
	# tar_load(sample_info)
	sample_info = add_blanks_pooled(colnames(cor_pca$correlation$cor), sample_info)
	cor_matrix = cor_pca$correlation$cor
	
	cor_order = data.frame(sample_id = colnames(cor_matrix))
	sample_info = dplyr::left_join(cor_order, sample_info, by = "sample_id")
	rownames(sample_info) = sample_info$sample_id
	
	treatment_order = visualizationQualityControl::similarity_reorderbyclass(cor_matrix, sample_classes = sample_info[, c("treatment"), drop = FALSE], transform = "sub_1")
	
	cor_matrix_treatment = cor_matrix[treatment_order$indices, treatment_order$indices]
	sample_info_treatment = sample_info[treatment_order$indices, ]
	
	median_cor_treatment = visualizationQualityControl::median_correlations(cor_matrix_treatment, sample_classes = sample_info_treatment$treatment)
	median_outlier_treatment = visualizationQualityControl::determine_outliers(median_correlations = median_cor_treatment)
	sample_info_treatment = dplyr::left_join(sample_info_treatment, median_outlier_treatment, by = "sample_id")
	
	low_range = round(min(cor_matrix_treatment), digits = 2)
	hi_range = round(max(cor_matrix_treatment), digits = 2)
	range_map = c(low_range, hi_range)
	n_color = 20
	colormap = circlize::colorRamp2(seq(range_map[1], range_map[2], length.out = n_color), viridis::viridis(n_color))
	
	annotate_treatment_groups = c("treatment", "outlier")
	
	row_annotation_treatment = ComplexHeatmap::HeatmapAnnotation(df = sample_info_treatment[, annotate_treatment_groups, drop = FALSE],
																												 col = color_scales, which = "row",
																												 show_annotation_name = FALSE)
	col_annotation_treatment = ComplexHeatmap::HeatmapAnnotation(df = sample_info_treatment[, annotate_treatment_groups, drop = FALSE],
																												 col = color_scales, which = "column",
																												 show_legend = FALSE, show_annotation_name = FALSE)
	
	treatment_heatmap = ComplexHeatmap::Heatmap(cor_matrix_treatment, col = colormap, name = "ICI-Kt",
																				bottom_annotation = col_annotation_treatment,
																				right_annotation = row_annotation_treatment,
																				cluster_rows = FALSE,
																				cluster_columns = FALSE,
																				column_title = "ICI-Kendell-tau Correlation Heatmap",
																				column_title_gp = grid::gpar(fontsize = 10),
																				column_names_gp = grid::gpar(fontsize = 8),
																				row_names_gp = grid::gpar(fontsize = 8))
	
	all_order = visualizationQualityControl::similarity_reorder(cor_matrix, transform = "sub_1")
	cor_matrix_all = cor_matrix[all_order$indices, all_order$indices]
	sample_info_all = sample_info[all_order$indices, ]
	
	median_cor_all = visualizationQualityControl::median_correlations(cor_matrix_all)
	median_outlier_all = visualizationQualityControl::determine_outliers(median_correlations = median_cor_all)
	sample_info_all = dplyr::left_join(sample_info_all, median_outlier_all, by = "sample_id")
	
	annotate_treatment_groups = c("treatment", "outlier")
	
	row_annotation_all = ComplexHeatmap::HeatmapAnnotation(df = sample_info_all[, annotate_treatment_groups, drop = FALSE],
																															 col = color_scales, which = "row",
																															 show_annotation_name = FALSE)
	col_annotation_all = ComplexHeatmap::HeatmapAnnotation(df = sample_info_all[, annotate_treatment_groups, drop = FALSE],
																															 col = color_scales, which = "column",
																															 show_legend = FALSE, show_annotation_name = FALSE)
	
	all_heatmap = ComplexHeatmap::Heatmap(cor_matrix_all, col = colormap, name = "ICI-Kt",
																							bottom_annotation = col_annotation_all,
																							right_annotation = row_annotation_all,
																							cluster_rows = FALSE,
																							cluster_columns = FALSE,
																							column_title = "ICI-Kendell-tau Correlation Heatmap",
																							column_title_gp = grid::gpar(fontsize = 10),
																							column_names_gp = grid::gpar(fontsize = 8),
																							row_names_gp = grid::gpar(fontsize = 8))
	
	tmp_all_pca = as.data.frame(cor_pca$all_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$all_pca$x))
	tmp_all_pca_var = visualizationQualityControl::visqc_score_contributions(cor_pca$all_pca$x)
	treatment_all_pca = dplyr::left_join(tmp_all_pca, sample_info_treatment, by = "sample_id")
	
	all_pca_plot = treatment_all_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(x = tmp_all_pca_var$labels[1], y = tmp_all_pca_var$labels[2])
	
	tmp_noblanks_pca = as.data.frame(cor_pca$noblanks_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$noblanks_pca$x))
	tmp_noblanks_pca_var = visualizationQualityControl::visqc_score_contributions(cor_pca$noblanks_pca$x)
	treatment_noblanks_pca = dplyr::left_join(tmp_noblanks_pca, sample_info_treatment, by = "sample_id")
	
	noblanks_pca_plot = treatment_noblanks_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(x = tmp_noblanks_pca_var$labels[1], y = tmp_noblanks_pca_var$labels[2])
	
	tmp_nooutlier_pca = as.data.frame(cor_pca$nooutlier_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$nooutlier_pca$x))
	tmp_nooutlier_pca_var = visualizationQualityControl::visqc_score_contributions(cor_pca$nooutlier_pca$x)
	treatment_nooutlier_pca = dplyr::left_join(tmp_nooutlier_pca, sample_info_treatment, by = "sample_id")
	
	nooutlier_pca_plot = treatment_nooutlier_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(x = tmp_nooutlier_pca_var$labels[1], y = tmp_nooutlier_pca_var$labels[2])
	
	median_cor_treatment_plot = sample_info_treatment |>
		ggplot(aes(x = treatment, y = med_cor, color = treatment, shape = outlier, group = treatment)) +
		geom_sina() +
		scale_color_manual(values = color_scales$treatment)
	
	list(heatmap = treatment_heatmap,
			 heatmap_all = all_heatmap,
			 median_cor = median_cor_treatment_plot,
			 pca_all = all_pca_plot,
			 pca_noblanks = noblanks_pca_plot,
			 pca_nooutlier = nooutlier_pca_plot)
		
}

sample_colors = function()
{
	treatment_colors = c(scale_color_discrete()$palette(4), "#000000")
	names(treatment_colors) = c("blank", "pooled", "normal_adjacent", "cancerous", "none")
	outlier_colors = scale_color_discrete()$palette(2)
	names(outlier_colors) = c("TRUE", "FALSE")
	list(treatment = treatment_colors,
			 outlier = outlier_colors)
}