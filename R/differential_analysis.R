determine_outliers = function(normalized_se) {
	# normalized_se = tar_read(lipidomics_keep)
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

	sample_info = colData(normalized_se)

	matrix_cor = sample_sample_cor$cor
	median_cor = visualizationQualityControl::median_correlations(
		matrix_cor,
		sample_info$treatment
	)
	outlier_frac = visualizationQualityControl::outlier_fraction(
		log1p(wide_matrix),
		sample_info$treatment,
		remove_missing = c(0, NA)
	)

	median_outliers = visualizationQualityControl::determine_outliers(
		median_correlations = median_cor,
		outlier_fraction = outlier_frac
	)

	rownames(median_outliers) = median_outliers$sample_id
	median_outliers = median_outliers[sample_info$sample_id, ]

	normalized_se[["outlier"]] = median_outliers$outlier
	normalized_se[["Outlier"]] = median_outliers$outlier
	normalized_se[["med_cor"]] = median_outliers$med_cor
	normalized_se[["frac"]] = median_outliers$frac
	normalized_se[["score"]] = median_outliers$score

	return(normalized_se)
}

collapse_deseq_replicates = function(outlier_se) {
	# outlier_se = tar_read(pm_outliers)
	outliers = outlier_se$outlier
	outlier_se = outlier_se[, !outliers]
	keep_reps = !is.na(outlier_se$replicate)
	outlier_se = outlier_se[, keep_reps]
	collapsed_se = collapseReplicates(
		outlier_se,
		outlier_se$replicate,
		outlier_se$sample_id,
		renameCols = TRUE
	)

	collapsed_se$treatment = forcats::fct_drop(collapsed_se$treatment)
	collapsed_se$patient = factor(collapsed_se$patient)
	collapsed_se$sample_id = colnames(counts(collapsed_se))
	collapsed_se
}

collapse_metabolomics_replicates = function(outlier_se) {
	# outlier_se = tar_read(bioamines_outliers)
	outliers = outlier_se$outlier
	pool_blanks = grepl("^pool|^blank", outlier_se$sample_id)
	outlier_se = outlier_se[, !(outliers | pool_blanks)]

	grouping_var = outlier_se$replicate

	grouping_list = split(colnames(outlier_se), grouping_var)
	n_in_each_group = purrr::map_int(grouping_list, length)

	if (max(n_in_each_group) == 1) {
		return(outlier_se)
	} else {
		org_counts = assays(outlier_se)$counts
		new_counts = purrr::imap(grouping_list, \(use_cols, new_id) {
			tmp_counts = org_counts[, use_cols, drop = FALSE]
			if (ncol(tmp_counts) == 1) {
				colnames(tmp_counts) = new_id
				return(tmp_counts)
			} else {
				out_counts = matrix(
					rowMeans(tmp_counts, na.rm = TRUE),
					nrow = nrow(tmp_counts),
					ncol = 1,
					byrow = TRUE
				)
				rownames(out_counts) = rownames(tmp_counts)
				colnames(out_counts) = new_id
				return(out_counts)
			}
		})
		new_counts = do.call(cbind, new_counts)

		new_coldata = colData(outlier_se)[colnames(new_counts), ]

		out_se = SummarizedExperiment(
			assays = list(counts = new_counts),
			rowData = rowData(outlier_se),
			colData = new_coldata
		)
		return(out_se)
	}
}


calculate_deseq_stats = function(
	rna_se,
	which = "treatment",
	contrast = c("treatment", "cancerous", "normal_adjacent"),
	fit_type = "local",
	named_only = FALSE
) {
	# rna_se = tar_read(rna_collapsed)
	# which = "treatment"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# rna_se = tar_read(rna_paired)
	# which = "patient"
	# fit_type = "parametric"
	#
	# rna_se = tar_read(pm_collapsed)
	# which = "treatment"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	#
	# rna_se = tar_read(lipidomics_collapsed)
	# which = "treatment"
	if (which %in% "treatment") {
		design(rna_se) = ~treatment
	} else {
		design(rna_se) = ~ patient + treatment
	}

	rna_deseq = DESeq(rna_se, fitType = fit_type)
	rna_results = results(rna_deseq, contrast = contrast) |> as.data.frame()
	rna_info = rowData(rna_se) |> as.data.frame()
	rna_results = cbind(rna_results, rna_info)

	if (named_only) {
		if ("metabolite_id" %in% colnames(rna_results)) {
			has_id = !is.na(rna_results[["metabolite_id"]])
			rna_results = rna_results[has_id, ]
			rna_results$padj = p.adjust(rna_results$pvalue)
		}
	}
	rna_results
}


filter_to_pairs = function(data_se) {
	# data_se = tar_read(rna_collapsed)
	patient_data = colData(data_se) |>
		as.data.frame() |>
		dplyr::select(sample_id, patient)
	n_patient_2 = patient_data |>
		dplyr::group_by(patient) |>
		dplyr::summarise(n_sample = dplyr::n()) |>
		dplyr::filter(n_sample == 2)
	keep_samples = data_se$patient %in% n_patient_2$patient
	data_se = data_se[, keep_samples]

	# remove the unused patients from the factor levels
	data_se$patient = factor(data_se$patient)

	sample_order = sample(length(data_se$patient), length(data_se$patient))
	data_se = data_se[, sample_order]

	return(data_se)
}

calculate_metabolomics_stats = function(
	data_se,
	paired = NULL,
	contrast = c("treatment", "cancerous", "normal_adjacent"),
	missing = c(0, NA)
) {
	# data_se = tar_read(pm_collapsed)
	# paired = NULL
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)

	# data_se = tar_read(pm_paired)
	# paired = "patient"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)

	# data_se = tar_read(bioamines_paired)
	# data_se = tar_read(bioamines_collapsed)
	# paired = "treatment"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)
	contrast_var = contrast[1]
	numerator = contrast[2]
	denominator = contrast[3]

	sample_info = colData(data_se) |> as.data.frame()
	# actually arrange things so they are paired
	t_paired = FALSE
	if (!is.null(paired) && (paired %in% "patient")) {
		sample_info = sample_info |>
			dplyr::arrange(patient, treatment)
		t_paired = TRUE
	}
	num_samples = sample_info[["sample_id"]][
		sample_info[[contrast_var]] == numerator
	]
	denom_samples = sample_info[["sample_id"]][
		sample_info[[contrast_var]] == denominator
	]

	counts = assays(data_se)$counts
	counts[counts %in% missing] = NA
	log_counts = log2(counts)
	use_missing = log2(sample_info$minimum)
	names(use_missing) = sample_info$sample_id
	num_counts = log_counts[, num_samples]
	denom_counts = log_counts[, denom_samples]

	out_stats = purrr::map(rownames(log_counts), \(in_row) {
		num_tmp = num_counts[in_row, ]
		denom_tmp = denom_counts[in_row, ]

		num_na = is.na(num_tmp)
		denom_na = is.na(denom_tmp)

		num_n = sum(!num_na)
		na_num_locs = names(num_tmp)[num_na]
		num_tmp[na_num_locs] = use_missing[na_num_locs]
		denom_n = sum(!denom_na)
		na_den_locs = names(denom_tmp)[denom_na]
		denom_tmp[na_den_locs] = use_missing[na_den_locs]
		if ((denom_n < 3) && (num_n < 3)) {
			return(NULL)
		}

		t_res = broom::tidy(t.test(num_tmp, denom_tmp, paired = t_paired))
		if (t_paired) {
			names(t_res)[1] = contrast[1]
		} else {
			names(t_res)[c(1, 2, 3)] = contrast
		}
		t_res[[paste0("n_", contrast[2])]] = num_n
		t_res[[paste0("n_", contrast[3])]] = denom_n
		t_res[["feature_id"]] = in_row
		return(t_res)
	}) |>
		purrr::list_rbind()
	out_stats$padj = p.adjust(out_stats$p.value, method = "BH")
	feature_info = rowData(data_se)[out_stats$feature_id, ] |> as.data.frame()

	out_stats = dplyr::left_join(out_stats, feature_info, by = "feature_id")
	out_stats$log2FoldChange = out_stats[[contrast[1]]]

	return(out_stats)
}

calculate_metabolomics_stats_aov = function(
	data_se,
	paired = NULL,
	contrast = c("treatment", "cancerous", "normal_adjacent"),
	missing = c(0, NA)
) {
	# data_se = tar_read(pr_collapsed)
	# paired = NULL
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)
	#
	# data_se = tar_read(pr_paired)
	# paired = "patient"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)

	# data_se = tar_read(bioamines_paired)
	# data_se = tar_read(bioamines_collapsed)
	# paired = NULL
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# missing = c(0, NA)
	contrast_var = contrast[1]
	numerator = contrast[2]
	denominator = contrast[3]

	sample_info = colData(data_se) |> as.data.frame()

	counts = assays(data_se)$counts
	counts[counts %in% missing] = NA
	log_counts = log2(counts)
	use_missing = signif(min(log_counts, na.rm = TRUE), digits = 2)

	use_df = sample_info[, c("treatment", "patient")]
	num_locs = use_df$treatment == numerator
	denom_locs = use_df$treatment == denominator

	out_stats = purrr::map(rownames(log_counts), \(in_row) {
		use_vals = log_counts[in_row, ]

		aov_setup = cbind(data.frame(abundance = use_vals), use_df)

		num_n = sum(!is.na(aov_setup$abundance[num_locs]))
		denom_n = sum(!is.na(aov_setup$abundance[denom_locs]))
		if (num_n < 3) {
			aov_setup$abundance[num_locs] = use_missing
		} else {
			num_na = is.na(aov_setup$abundance) & num_locs
			aov_setup$abundance[num_na] = use_missing
		}
		if (denom_n < 3) {
			aov_setup$abundance[denom_locs] = use_missing
		} else {
			denom_na = is.na(aov_setup$abundance) & denom_locs
			aov_setup$abundance[denom_na] = use_missing
		}

		if ((denom_n < 3) && (num_n < 3)) {
			return(NULL)
		}

		if (is.null(paired) || (paired == "paired")) {
			aov_res = broom::tidy(aov(
				abundance ~ treatment,
				data = aov_setup
			)) |>
				dplyr::filter(term %in% "treatment")
		} else {
			aov_res = broom::tidy(aov(
				abundance ~ patient + treatment,
				data = aov_setup
			)) |>
				dplyr::filter(term %in% "treatment")
		}
		diff_vals = data.frame(
			numerator = mean(aov_setup$abundance[
				aov_setup$treatment %in% numerator
			]),
			denominator = mean(aov_setup$abundance[
				aov_setup$treatment %in% denominator
			])
		)
		diff_vals$contrast = diff_vals$numerator - diff_vals$denominator
		names(diff_vals) = c(numerator, denominator, contrast_var)
		out_res = cbind(aov_res, diff_vals)
		out_res$feature_id = in_row
		out_res
	}) |>
		purrr::list_rbind()

	out_stats$padj = p.adjust(out_stats$p.value, method = "BH")
	feature_info = rowData(data_se)[out_stats$feature_id, ] |> as.data.frame()

	out_stats = dplyr::left_join(out_stats, feature_info, by = "feature_id")
	out_stats$log2FoldChange = out_stats[[contrast[1]]]

	return(out_stats)
}

calculate_limma_stats = function(
	data_se,
	which = "patient",
	contrast = c("treatment", "normal_adjacent", "cancerous")
) {
	# data_se = tar_read(lipidomics_paired)
	# which = "patient"
	# contrast = c("treatment", "normal_adjacent",
	# 						"cancerous")

	sample_info = colData(data_se) |> tibble::as_tibble()
	sample_info$treatment = factor(
		sample_info$treatment,
		levels = contrast[c(2, 3)]
	)
	design_matrix = model.matrix(~ patient + treatment, data = sample_info)
	ncol_design = ncol(design_matrix)
	colnames(design_matrix)[ncol_design] = "treatment"

	gene_data = counts(data_se, normalized = FALSE)
	gene_data[gene_data == 0] = NA
	sample_medians = apply(gene_data, 2, median, na.rm = TRUE)
	matrix_medians = matrix(
		sample_medians,
		nrow = nrow(gene_data),
		ncol = ncol(gene_data),
		byrow = TRUE
	)

	gene_norm = gene_data / matrix_medians

	log2_norm = log2(gene_norm)

	lm_fit = limma::lmFit(log2_norm, design_matrix)
	e_fit = limma::eBayes(lm_fit)

	results = limma::topTable(
		e_fit,
		coef = "treatment",
		number = Inf,
		p.value = 1
	)
	results$feature_id = rownames(results)
	results
}

merge_list = function(in_data) {
	merged_vals = purrr::imap(in_data, \(x, id) {
		x |>
			dplyr::mutate(type = id)
	}) |>
		purrr::list_rbind()
	merged_vals
}

map_metabolomics_chebi = function(in_data, chebi_inchikey) {
	# in_data = tar_read(metabolomics_de_treatment_list)
	# tar_load(chebi_inchikey)

	de_vals = dplyr::left_join(
		in_data,
		chebi_inchikey,
		by = c("in_chi_key" = "in_ch_i_key"),
		relationship = "many-to-many"
	)
	de_vals = de_vals |>
		dplyr::mutate(
			feature_org = feature_id,
			feature_id = as.character(chebi_id)
		)
	de_vals
}

map_metabolomics_kegg = function(in_data, inchikey_kegg) {
	# in_data = tar_read(metabolomics_de_treatment_list)
	# tar_load(inchikey_kegg)

	de_vals = dplyr::left_join(
		in_data |> dplyr::select(-kegg),
		inchikey_kegg,
		by = "in_chi_key",
		relationship = "many-to-many"
	)
	de_vals = de_vals |>
		dplyr::mutate(
			feature_org = feature_id,
			feature_id = as.character(kegg)
		) |>
		dplyr::filter(!is.na(feature_id))
	de_vals
}

compare_treatment_patient = function(de_treatment, de_patient) {
	# de_treatment = tar_read(metabolomics_de_treatment)
	# de_patient = tar_read(metabolomics_de_patient)
	#
	# de_treatment = tar_read(rna_de_treatment)
	# de_patient = tar_read(rna_de_patient)
	#
	# de_treatment = tar_read(metabolomics_de_aov_treatment)
	# de_patient = tar_read(metabolomics_de_aov_patient)
	#

	de_treatment_mod = de_treatment |>
		dplyr::select(log2FoldChange, pvalue, padj, feature_id) |>
		dplyr::distinct()
	de_patient_mod = de_patient |>
		dplyr::select(log2FoldChange, pvalue, padj, feature_id) |>
		dplyr::distinct()
	de_merge = dplyr::left_join(
		de_treatment_mod,
		de_patient_mod,
		by = "feature_id",
		suffix = c(".treatment", ".patient")
	)

	de_merge
}


information_volume <- function(in_x, in_y) {
	in_both <- sum(in_x & in_y)
	in_both
}

both_samples <- function(in_x, in_y) {
	in_both <- sum(in_x & in_y)
	in_both
}

either_samples <- function(in_x, in_y) {
	in_either <- sum(in_x | in_y)
	in_either
}

neither_sample <- function(in_x, in_y) {
	not_in_both <- sum(!in_x & !in_y)
	not_in_both
}

calculate_jaccard <- function(in_both, in_either) {
	in_both / in_either
}

calculate_information_consistency <- function(
	in_both,
	not_in_both,
	total_features
) {
	(in_both + not_in_both) / total_features
}

information_stats = function(num_na, denom_na) {}

extract_deseq_patient_logratios = function(
	rna_se,
	contrast = c("treatment", "cancerous", "normal_adjacent"),
	data_type = "rna"
) {
	# rna_se = tar_read(rna_paired)
	# contrast = c("treatment", "cancerous",
	# 						 "normal_adjacent")

	norm_counts = DESeq2::counts(rna_se, normalized = TRUE)

	# min_value = min(norm_counts[norm_counts > 0])
	# norm_counts[norm_counts == 0] = min_value

	log_counts = log2(norm_counts)
	sample_info = colData(rna_se) |> tibble::as_tibble()

	samples_by_patient = split(sample_info, sample_info$patient)

	log_ratios = purrr::map(samples_by_patient, \(in_patient) {
		# in_patient = samples_by_patient[[5]]
		use_column = contrast[1]
		ref_sample = in_patient |>
			dplyr::filter(treatment %in% contrast[3]) |>
			dplyr::pull(replicate)
		treat_sample = in_patient |>
			dplyr::filter(treatment %in% contrast[2]) |>
			dplyr::pull(replicate)
		log_counts[, treat_sample] - log_counts[, ref_sample]
	}) |>
		dplyr::bind_cols() |>
		as.matrix()
	log_ratios[
		is.infinite(log_ratios) | is.na(log_ratios) | is.nan(log_ratios)
	] = NA
	rownames(log_ratios) = rownames(log_counts)

	log_ratios
}


create_logratio_heatmap = function(
	rna_patient_logratios,
	rna_de_patient,
	max_padj = 0.01,
	label = "Genes"
) {
	# tar_load(c(rna_patient_logratios,
	# 					 rna_de_patient))
	# max_padj = 0.01

	force(max_padj)
	rna_sig = rna_de_patient |>
		dplyr::filter(padj <= max_padj)

	rna_sig_logratio = rna_patient_logratios[rna_sig$feature_id, ]

	n_miss_gene = rowSums(is.na(rna_sig_logratio))
	rna_sig_logratio = rna_sig_logratio[n_miss_gene <= 2, ]
	use_range = range(rna_sig_logratio, na.rm = TRUE)
	use_range[1] = ceiling(use_range[1])
	use_range[2] = floor(use_range[2])
	n_value = 20
	heatmap_colors = circlize::colorRamp2(
		seq(use_range[1], use_range[2], length.out = n_value),
		scico::scico(n_value, palette = "berlin")
	)

	out_heatmap = Heatmap(
		rna_sig_logratio,
		col = heatmap_colors,
		name = "Log2FC",
		show_row_names = FALSE,
		show_column_names = FALSE,
		row_title = "Features",
		column_title = "Patients",
		row_labels = fix_labels(rownames(rna_sig_logratio))
	)

	out_heatmap
}

check_gene_rankings = function(
	rna_de_patient,
	rna_patient_logratios,
	features = c("HBA", "BASP1", "OLFML3", "HK3", "AKR1B10", "PEG10", "MYBL2"),
	n_missing = 0
) {
	# tar_load(c(rna_de_patient, rna_patient_logratios))
	# features = c("HBA", "BASP1", "OLFML3", "HK3", "AKR1B10", "PEG10", "MYBL2")
	# n_missing = 0
	max_padj = 0.01
	n_miss_logratios = data.frame(
		feature_id = rownames(rna_patient_logratios),
		n_na = rowSums(is.na(rna_patient_logratios))
	)
	n_miss_low = n_miss_logratios |>
		dplyr::filter(n_na <= n_missing)
	keep_patient = rna_de_patient |>
		dplyr::filter(
			padj <= max_padj,
			feature_id %in% n_miss_low$feature_id,
			biotype %in% "protein_coding"
		)

	pos_ranking = keep_patient |>
		dplyr::filter(log2FoldChange > 0) |>
		dplyr::arrange((log2FoldChange)) |>
		dplyr::mutate(rank = rank(dplyr::desc(log2FoldChange)))
	neg_ranking = keep_patient |>
		dplyr::filter(log2FoldChange < 0) |>
		dplyr::arrange(log2FoldChange) |>
		dplyr::mutate(rank = rank(log2FoldChange))

	all_ranking = dplyr::bind_rows(pos_ranking, neg_ranking)

	all_ranking = all_ranking |>
		dplyr::mutate(
			is_interesting = dplyr::case_when(
				name %in% features ~ "yes",
				TRUE ~ "no"
			)
		)
	logratio_df = rna_patient_logratios |>
		as.data.frame() |>
		tibble::rownames_to_column(var = "feature_id")

	all_ranking = dplyr::left_join(all_ranking, logratio_df, by = "feature_id")
	write.table(
		all_ranking,
		file = "docs/helsley_zeromissing_ranks.csv",
		sep = ",",
		row.names = FALSE,
		col.names = TRUE
	)
}

create_logratio_specific_heatmap = function(
	csv_file = "data/Top15UpandDown_Helsley.csv",
	fontsize = 10,
	label = "Genes"
) {
	# csv_file = "data/Top15UpandDown_Helsley.csv"
	# fontsize = 12
	# label = "Genes"
	csv_data = readr::read_csv(csv_file)

	csv_specific = csv_data |>
		dplyr::filter(is_interesting %in% "yes")

	rank_logratio = csv_specific |>
		dplyr::arrange(dplyr::desc(log2FoldChange))

	all_lfc = rank_logratio |>
		dplyr::select(tidyselect::starts_with("patient")) |>
		as.matrix()
	rownames(all_lfc) = rank_logratio$name

	max_range = ceiling(max(abs(all_lfc), na.rm = TRUE))
	use_range = c(-1 * max_range, max_range)
	n_value = 20
	heatmap_colors = circlize::colorRamp2(
		seq(use_range[1], use_range[2], length.out = n_value),
		scico::scico(n_value, palette = "berlin")
	)

	out_heatmap = Heatmap(
		all_lfc,
		col = heatmap_colors,
		name = "Log2FC",
		show_row_names = TRUE,
		show_column_names = FALSE,
		row_title = label,
		column_title = "Patients",
		cluster_rows = FALSE,
		row_names_gp = gpar(fontsize = fontsize),
		row_labels = fix_labels(rownames(all_lfc))
	)
	return(out_heatmap)
}

create_logratio_heatmap_small = function(
	rna_patient_logratios,
	rna_de_patient,
	limit = 15,
	n_missing = 2,
	id = "name",
	fontsize = 10,
	label = "Genes"
) {
	# tar_load(c(rna_patient_logratios,
	# 					 rna_de_patient))
	# limit = 15
	# id = "name"
	#
	# rna_patient_logratios = tar_read(lipidomics_patient_logratios)
	# rna_de_patient = tar_read(lipidomics_de_patient)
	# limit = 15
	# id = "metabolite_id"

	# tar_load(c(rna_patient_logratios,
	# 					 rna_de_patient))
	# limit = 15
	# id = "name"
	# n_missing = 0

	max_padj = 0.01
	n_miss_logratios = data.frame(
		feature_id = rownames(rna_patient_logratios),
		n_na = rowSums(is.na(rna_patient_logratios))
	)
	n_miss_low = n_miss_logratios |>
		dplyr::filter(n_na <= n_missing)
	keep_patient = rna_de_patient |>
		dplyr::filter(padj <= max_padj, feature_id %in% n_miss_low$feature_id)

	pos_15 = keep_patient |>
		dplyr::filter(log2FoldChange > 0) |>
		dplyr::arrange(dplyr::desc(log2FoldChange)) |>
		dplyr::slice_head(n = limit)
	neg_15 = keep_patient |>
		dplyr::filter(log2FoldChange < 0) |>
		dplyr::arrange(log2FoldChange) |>
		dplyr::slice_head(n = limit) |>
		dplyr::arrange(dplyr::desc(log2FoldChange))

	all_15 = dplyr::bind_rows(pos_15, neg_15)

	all_lfc = rna_patient_logratios[all_15$feature_id, ]

	max_range = ceiling(max(abs(all_lfc), na.rm = TRUE))
	use_range = c(-1 * max_range, max_range)
	n_value = 20
	heatmap_colors = circlize::colorRamp2(
		seq(use_range[1], use_range[2], length.out = n_value),
		scico::scico(n_value, palette = "berlin")
	)

	out_heatmap = Heatmap(
		all_lfc,
		col = heatmap_colors,
		name = "Log2FC",
		show_row_names = TRUE,
		show_column_names = FALSE,
		row_title = label,
		column_title = "Patients",
		cluster_rows = FALSE,
		row_labels = fix_labels(all_15[[id]]),
		row_names_gp = gpar(fontsize = fontsize)
	)
	return(out_heatmap)
}
