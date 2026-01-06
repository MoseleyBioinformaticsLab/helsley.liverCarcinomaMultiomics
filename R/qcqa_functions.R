rename_experimental_samples = function(metabolomics_dataset) {
	# metabolomics_dataset = tar_read(bioamines)
	names(metabolomics_dataset) = gsub("^x", "s", names(metabolomics_dataset))
	metabolomics_dataset
}

median_normalization = function(feature_se, min_type = "global") {
	# feature_se = tar_read(bioamines)
	# feature_se = tar_read(lipidomics)
	# feature_se = tar_read(primary_metabolism)
	feature_abundance = assays(feature_se)$counts
	sample_medians = matrixStats::colMedians(feature_abundance, na.rm = TRUE)
	median_matrix = matrix(
		sample_medians,
		nrow = nrow(feature_abundance),
		ncol = ncol(feature_abundance),
		byrow = TRUE
	)
	norm_abundance = feature_abundance / median_matrix

	assays(feature_se)$counts = norm_abundance
	feature_se$normalization = sample_medians

	if (min_type %in% "global") {
		min_vals = rep(min(norm_abundance, na.rm = TRUE), length(sample_medians))
		names(min_vals) = names(sample_medians)
	} else if (min_type %in% "sample") {
		min_all = min(feature_abundance, na.rm = TRUE)
		min_vals = min_all * sample_medians
	} else if (min_type %in% "patient") {
		split_sample_patient = split(feature_se$sample_id, feature_se$patient)
		names(split_sample_patient) = NULL
		min_patient = purrr::map(split_sample_patient, \(in_patients) {
			min_value = min(norm_abundance[, in_patients], na.rm = TRUE)
			min_out = rep(min_value, length(in_patients))
			names(min_out) = in_patients
			min_out
		})
		min_vals = unlist(min_patient)
		min_vals = min_vals[names(sample_medians)]
	}

	feature_se$minimum = min_vals

	feature_se
}

matrix_2_long = function(wide_matrix) {
	# wide_matrix = tar_load(rna_norm)
	tmp_df = tibble::as_tibble(wide_matrix) |>
		dplyr::mutate(feature_id = rownames(wide_matrix))
	long_df = tmp_df |>
		tidyr::pivot_longer(
			cols = tidyselect::where(is.numeric),
			names_to = "sample_id",
			values_to = "abundance"
		)
	long_df
}

long_2_matrix = function(long_data) {
	# long_data = tar_read(rna_norm)
	wide_df = long_data |>
		tidyr::pivot_wider(
			id_cols = feature_id,
			names_from = "sample_id",
			values_from = "abundance"
		)
	wide_matrix = as.matrix(wide_df |> dplyr::select(-feature_id))
	rownames(wide_matrix) = wide_df$feature_id
	wide_matrix
}

keep_presence = function(normalized_se, fraction = 0.75) {
	# normalized_se = tar_read(bioamines_norm)
	# fraction = 0.75

	keep_norm = visualizationQualityControl::keep_non_missing_percentage(
		assays(normalized_se)$counts,
		colData(normalized_se)$treatment,
		keep_num = fraction,
		missing_value = c(NA, 0)
	)
	out_se = normalized_se[keep_norm, ]
	return(out_se)
}

floor_values = function(rna_matrix) {
	new_matrix = rna_matrix
	new_matrix[rna_matrix < 10] = 0
	new_matrix
}


keep_presence_dds = function(rna_dds, fraction = 0.75) {
	# rna_dds = tar_read(rna_dds)
	# fraction = 0.75
	norm_wide = DESeq2::counts(rna_dds)
	sample_info = as.data.frame(colData(rna_dds))
	# set anything less than 10 to 0 so it still counts as missing
	norm_wide[norm_wide < 10] = 0
	keep_norm = visualizationQualityControl::keep_non_missing_percentage(
		norm_wide,
		sample_info$treatment,
		keep_num = fraction,
		missing_value = c(NA, 0)
	)
	rna_dds = rna_dds[keep_norm, ]
	return(rna_dds)
}

sample_correlations_pca = function(normalized_se) {
	# normalized_se = tar_read(bioamines_keep)
	# normalized_se = tar_read(rna_keep)
	# normalized_se = tar_read(primary_metabolism_keep)
	if (inherits(normalized_se, "DESeqDataSet")) {
		wide_matrix = DESeq2::counts(normalized_se, normalized = TRUE)
	} else {
		wide_matrix = assays(normalized_se)$counts
	}
	wide_matrix[is.nan(wide_matrix) | is.infinite(wide_matrix)] = NA

	sample_sample_cor = suppressWarnings(ICIKendallTau::ici_kendalltau(
		wide_matrix,
		global_na = c(0, NA),
		perspective = "global",
		scale_max = TRUE,
		diag_good = TRUE
	))

	if (any(is.na(wide_matrix))) {
		wide_matrix[is.na(wide_matrix)] = 0
	}

	sample_all_pca = prcomp(t(log1p(wide_matrix)), center = TRUE, scale. = FALSE)
	pool_blank_samples = grepl("^pool|^blank", colnames(wide_matrix))
	sample_noblanks_pca = prcomp(
		t(log1p(wide_matrix[, !pool_blank_samples])),
		center = TRUE,
		scale. = FALSE
	)

	sample_info = colData(normalized_se) |> as.data.frame()

	matrix_cor = sample_sample_cor$cor
	median_cor = visualizationQualityControl::median_correlations(
		matrix_cor,
		sample_info$treatment
	)
	outlier_frac = visualizationQualityControl::outlier_fraction(
		log1p(matrix_cor),
		sample_info$treatment,
		remove_missing = c(0, NA)
	)
	median_outliers = visualizationQualityControl::determine_outliers(
		median_correlations = median_cor,
		outlier_fraction = outlier_frac
	)
	sample_info = dplyr::left_join(sample_info, median_outliers, by = "sample_id")
	sample_info = sample_info |>
		dplyr::mutate(Outlier = outlier)
	use_samples = median_outliers |>
		dplyr::filter(!outlier, !grepl("blank|pool", sample_id)) |>
		dplyr::pull(sample_id)

	sample_nooutlier_pca = prcomp(
		t(log1p(wide_matrix[, use_samples])),
		center = TRUE,
		scale. = FALSE
	)

	return(list(
		correlation = sample_sample_cor,
		all_pca = sample_all_pca,
		noblanks_pca = sample_noblanks_pca,
		nooutlier_pca = sample_nooutlier_pca,
		info = sample_info
	))
}

create_qcqa_plots = function(cor_pca, color_scales) {
	# cor_pca = tar_read(bioamines_cor_pca)
	# cor_pca = tar_read(primary_metabolism_cor_pca)
	# cor_pca = tar_read(rna_cor_pca)
	# tar_load(color_scales)
	sample_info = cor_pca$info
	cor_matrix = cor_pca$correlation$cor

	rownames(sample_info) = sample_info$sample_id

	treatment_order = visualizationQualityControl::similarity_reorderbyclass(
		cor_matrix,
		sample_classes = sample_info[, c("treatment"), drop = FALSE],
		transform = "sub_1"
	)

	cor_matrix_treatment = cor_matrix[
		treatment_order$indices,
		treatment_order$indices
	]
	sample_info_treatment = sample_info[treatment_order$indices, ]

	low_range = round(min(cor_matrix_treatment), digits = 2)
	hi_range = round(max(cor_matrix_treatment), digits = 2)
	range_map = c(low_range, hi_range)
	n_color = 20
	colormap = circlize::colorRamp2(
		seq(range_map[1], range_map[2], length.out = n_color),
		viridis::viridis(n_color)
	)

	annotate_treatment_groups = c("treatment", "outlier")

	row_annotation_treatment = ComplexHeatmap::HeatmapAnnotation(
		df = sample_info_treatment[, annotate_treatment_groups, drop = FALSE],
		col = color_scales,
		which = "row",
		show_annotation_name = FALSE
	)
	col_annotation_treatment = ComplexHeatmap::HeatmapAnnotation(
		df = sample_info_treatment[, annotate_treatment_groups, drop = FALSE],
		col = color_scales,
		which = "column",
		show_legend = FALSE,
		show_annotation_name = FALSE
	)

	treatment_heatmap = ComplexHeatmap::Heatmap(
		cor_matrix_treatment,
		col = colormap,
		name = "ICI-Kt",
		bottom_annotation = col_annotation_treatment,
		right_annotation = row_annotation_treatment,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		column_title = "ICI-Kendell-tau Correlation Heatmap",
		column_title_gp = grid::gpar(fontsize = 10),
		column_names_gp = grid::gpar(fontsize = 8),
		row_names_gp = grid::gpar(fontsize = 8)
	)

	all_order = visualizationQualityControl::similarity_reorder(
		cor_matrix,
		transform = "sub_1"
	)
	cor_matrix_all = cor_matrix[all_order$indices, all_order$indices]
	sample_info_all = sample_info[all_order$indices, ]

	median_cor_all = visualizationQualityControl::median_correlations(
		cor_matrix_all
	)
	median_outlier_all = visualizationQualityControl::determine_outliers(
		median_correlations = median_cor_all
	)
	tmp_all = c(
		"sample_id",
		base::setdiff(
			names(sample_info),
			base::intersect(names(sample_info), names(median_outlier_all))
		)
	)
	sample_info_all = dplyr::left_join(
		sample_info_all[, tmp_all],
		median_outlier_all,
		by = "sample_id"
	)

	annotate_treatment_groups = c("treatment", "outlier")

	row_annotation_all = ComplexHeatmap::HeatmapAnnotation(
		df = sample_info_all[, annotate_treatment_groups, drop = FALSE],
		col = color_scales,
		which = "row",
		show_annotation_name = FALSE
	)
	col_annotation_all = ComplexHeatmap::HeatmapAnnotation(
		df = sample_info_all[, annotate_treatment_groups, drop = FALSE],
		col = color_scales,
		which = "column",
		show_legend = FALSE,
		show_annotation_name = FALSE
	)

	all_heatmap = ComplexHeatmap::Heatmap(
		cor_matrix_all,
		col = colormap,
		name = "ICI-Kt",
		bottom_annotation = col_annotation_all,
		right_annotation = row_annotation_all,
		cluster_rows = FALSE,
		cluster_columns = FALSE,
		column_title = "ICI-Kendell-tau Correlation Heatmap",
		column_title_gp = grid::gpar(fontsize = 10),
		column_names_gp = grid::gpar(fontsize = 8),
		row_names_gp = grid::gpar(fontsize = 8)
	)

	tmp_all_pca = as.data.frame(cor_pca$all_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$all_pca$x))
	tmp_all_pca_var = visualizationQualityControl::visqc_score_contributions(
		cor_pca$all_pca$x
	)
	treatment_all_pca = dplyr::left_join(
		tmp_all_pca,
		sample_info_treatment,
		by = "sample_id"
	)

	all_pca_plot = treatment_all_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(x = tmp_all_pca_var$labels[1], y = tmp_all_pca_var$labels[2])

	tmp_noblanks_pca = as.data.frame(cor_pca$noblanks_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$noblanks_pca$x))
	tmp_noblanks_pca_var = visualizationQualityControl::visqc_score_contributions(
		cor_pca$noblanks_pca$x
	)
	treatment_noblanks_pca = dplyr::left_join(
		tmp_noblanks_pca,
		sample_info_treatment,
		by = "sample_id"
	)

	noblanks_pca_plot = treatment_noblanks_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(x = tmp_noblanks_pca_var$labels[1], y = tmp_noblanks_pca_var$labels[2])

	tmp_nooutlier_pca = as.data.frame(cor_pca$nooutlier_pca$x) |>
		dplyr::mutate(sample_id = rownames(cor_pca$nooutlier_pca$x))
	tmp_nooutlier_pca_var = visualizationQualityControl::visqc_score_contributions(
		cor_pca$nooutlier_pca$x
	)
	treatment_nooutlier_pca = dplyr::left_join(
		tmp_nooutlier_pca,
		sample_info_treatment,
		by = "sample_id"
	)

	tmp_nooutlier_pca_anova = visualizationQualityControl::visqc_test_pca_scores(
		as.matrix(treatment_nooutlier_pca[, 1:15]),
		treatment_nooutlier_pca[, c("treatment", "patient")]
	)

	nooutlier_pca_plot = treatment_nooutlier_pca |>
		ggplot(aes(x = PC1, y = PC2, color = treatment, shape = outlier)) +
		geom_point(size = 2) +
		scale_color_manual(values = color_scales$treatment) +
		labs(
			x = tmp_nooutlier_pca_var$labels[1],
			y = tmp_nooutlier_pca_var$labels[2]
		)

	median_cor_treatment_plot = sample_info_treatment |>
		ggplot(aes(
			x = treatment,
			y = med_cor,
			color = treatment,
			shape = outlier,
			group = treatment
		)) +
		geom_sina() +
		scale_color_manual(values = color_scales$treatment)

	list(
		heatmap = treatment_heatmap,
		heatmap_all = all_heatmap,
		median_cor = sample_info_treatment,
		pca_all = all_pca_plot,
		pca_noblanks = noblanks_pca_plot,
		pca_nooutlier = nooutlier_pca_plot,
		pca_nooutlier_values = treatment_nooutlier_pca,
		pca_nooutlier_variance = tmp_nooutlier_pca_var,
		pca_nooutlier_anova = tmp_nooutlier_pca_anova
	)
}

create_color_scales = function() {
	treatment_colors = c(scale_color_discrete()$palette(4), "#000000")
	names(treatment_colors) = c(
		"blank",
		"pooled",
		"normal_adjacent",
		"cancerous",
		"none"
	)
	outlier_colors = scale_color_discrete()$palette(2)
	names(outlier_colors) = c("TRUE", "FALSE")

	just_normal_cancer = c("Normal Adjacent" = "#84B44C", "Cancerous" = "#CF7ABF")
	treatment_normal_cancer = c(
		"normal_adjacent" = "#84B44C",
		"cancerous" = "#CF7ABF"
	)
	volcano_plots = c("negative" = "#5DA5DD", "positive" = "#C6776C")

	list(
		treatment = treatment_colors,
		outlier = outlier_colors,
		normal_cancer = just_normal_cancer,
		treatment_normal_cancer = treatment_normal_cancer,
		volcano = volcano_plots
	)
}

get_n_features = function(in_data) {
	if (inherits(in_data, "data.frame")) {
		return(length(unique(in_data$feature_id)))
	} else {
		return(nrow(in_data))
	}
}

check_left_censoring = function(dds_obj) {
	left_test = ICIKendallTau::test_left_censorship(
		counts(dds_obj),
		sample_classes = dds_obj$treatment
	)

	left_ranks = ICIKendallTau::rank_order_data(
		counts(dds_obj),
		sample_classes = dds_obj$treatment
	)

	list(test = left_test, ranks = left_ranks)
}

count_n_samples = function(dds_obj, id) {
	# dds_obj = tar_read(rna_collapsed)
	# id = "Transcriptomics"
	n_treatment = colData(dds_obj) |>
		tibble::as_tibble() |>
		dplyr::filter(!outlier) |>
		dplyr::select(replicate, Treatment) |>
		dplyr::distinct() |>
		dplyr::group_by(Treatment) |>
		dplyr::summarise(n_sample = dplyr::n()) |>
		dplyr::mutate(Dataset = id) |>
		tidyr::pivot_wider(names_from = "Treatment", values_from = "n_sample")
	n_treatment
}

calculate_variances = function(dds_obj) {
	# dds_obj = tar_read(rna_paired)
	norm_counts = counts(dds_obj, normalized = TRUE, replaced = FALSE) |>
		as.data.frame() |>
		tibble::rownames_to_column("feature_id")
	long_counts = norm_counts |>
		tidyr::pivot_longer(
			-feature_id,
			names_to = "sample_id",
			values_to = "intensity"
		) |>
		dplyr::mutate(
			intensity = dplyr::case_when(
				intensity == 0 ~ NA,
				TRUE ~ intensity
			)
		)
	sample_info = colData(dds_obj) |> tibble::as_tibble()

	long_counts = dplyr::left_join(
		long_counts,
		sample_info[, c("sample_id", "treatment")],
		by = "sample_id"
	)
	variance_counts = long_counts |>
		dplyr::group_by(treatment, feature_id) |>
		dplyr::summarise(
			sd_intensity = sd(intensity, na.rm = TRUE),
			mean_intensity = mean(intensity, na.rm = TRUE),
			rsd_intensity = sd_intensity / mean_intensity,
			feature_id = feature_id[1]
		) |>
		dplyr::ungroup()

	var_ratios = variance_counts |>
		dplyr::group_by(feature_id) |>
		dplyr::arrange(treatment) |>
		dplyr::summarise(
			sd = log10(sd_intensity[2]) - log10(sd_intensity[1]),
			rsd = log10(rsd_intensity[2]) - log10(rsd_intensity[1])
		) |>
		dplyr::ungroup()

	long_var = var_ratios |>
		tidyr::pivot_longer(-feature_id, values_to = "ratio", names_to = "which")
	long_var
}
