determine_outliers = function(normalized_se)
{
	# normalized_se = tar_read(lipidomics_keep)
	# normalized_se = tar_read(rna_keep)
	# normalized_se = tar_read(primary_metabolism_keep)
	if (inherits(normalized_se, "DESeqDataSet")) {
		wide_matrix = DESeq2::counts(normalized_se, normalized = TRUE)
	} else {
		wide_matrix = assays(normalized_se)$counts
	}
	wide_matrix[is.nan(wide_matrix) | is.infinite(wide_matrix)] = NA
	
	sample_sample_cor = suppressWarnings(ICIKendallTau::ici_kendalltau(wide_matrix, global_na = c(0, NA), perspective = "global", scale_max = TRUE, diag_good = TRUE))
	
	if (any(is.na(wide_matrix))) {
		wide_matrix[is.na(wide_matrix)] = 0
	}
	
	sample_info = colData(normalized_se)
	
	matrix_cor = sample_sample_cor$cor
	median_cor = visualizationQualityControl::median_correlations(matrix_cor, sample_info$treatment)
	outlier_frac = visualizationQualityControl::outlier_fraction(log1p(wide_matrix), sample_info$treatment, remove_missing = c(0, NA))
	
	median_outliers = visualizationQualityControl::determine_outliers(median_correlations = median_cor, outlier_fraction = outlier_frac)
	
	rownames(median_outliers) = median_outliers$sample_id
	median_outliers = median_outliers[sample_info$sample_id, ]
	
	normalized_se[["outlier"]] = median_outliers$outlier
	normalized_se[["med_cor"]] = median_outliers$med_cor
	normalized_se[["frac"]] = median_outliers$frac
	normalized_se[["score"]] = median_outliers$score
	
	return(normalized_se)
}

collapse_deseq_replicates = function(outlier_se)
{
	# outlier_se = tar_read(rna_outliers)
	outliers = outlier_se$outlier
	outlier_se = outlier_se[, !outliers]
	collapsed_se = collapseReplicates(outlier_se, outlier_se$replicate, outlier_se$sample_id, renameCols = TRUE)
	collapsed_se
}

collapse_metabolomics_replicates = function(outlier_se)
{
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
		new_counts = purrr::imap(grouping_list, \(use_cols, new_id){
			tmp_counts = org_counts[, use_cols, drop = FALSE]
			if (ncol(tmp_counts) == 1) {
				colnames(tmp_counts) = new_id
				return(tmp_counts)
			}  else {
				out_counts = matrix(rowMeans(tmp_counts, na.rm = TRUE), nrow = nrow(tmp_counts), ncol = 1, byrow = TRUE)
				rownames(out_counts) = rownames(tmp_counts)
				colnames(out_counts) = new_id
				return(out_counts)
			}
		})
		new_counts = do.call(cbind, new_counts)
		
		new_coldata = colData(outlier_se)[colnames(new_counts), ]
		
		out_se = SummarizedExperiment(assays = list(counts = new_counts),
																	rowData = rowData(outlier_se),
																	colData = new_coldata)
		return(out_se)
		
	}
}


calculate_deseq_stats = function(rna_se,
																 which = "treatment",
																 contrast = c("treatment", "cancerous", "normal_adjacent"))
{
	# rna_se = tar_read(rna_collapsed)
	# which = "treatment"
	# contrast = c("treatment", "cancerous", "normal_adjacent")
	# rna_se = tar_read(rna_paired)
	# which = "patient"
	if (which %in% "treatment") {
		design(rna_se) = ~ treatment
	} else {
		design(rna_se) = ~ patient + treatment
		
	}
	
	rna_deseq = DESeq(rna_se)
	rna_results = results(rna_deseq, contrast = contrast) |> as.data.frame()
	rna_info = rowData(rna_se) |> as.data.frame()
	rna_results = cbind(rna_results, rna_info)
	rna_results
}

filter_to_pairs = function(data_se)
{
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

calculate_metabolomics_stats = function(data_se,
																				paired = NULL,
																				contrast = c("treatment", "cancerous", "normal_adjacent"),
																				missing = c(0, NA))
{
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
	num_samples = sample_info[["sample_id"]][sample_info[[contrast_var]] == numerator]
	denom_samples = sample_info[["sample_id"]][sample_info[[contrast_var]] == denominator]
	
	counts = assays(data_se)$counts
	counts[counts %in% missing] = NA
	log_counts = log2(counts)
	use_missing = signif(min(log_counts, na.rm = TRUE), digits = 2)
	num_counts = log_counts[, num_samples]
	denom_counts = log_counts[, denom_samples]
		
	out_stats = purrr::map(rownames(log_counts), \(in_row){
		num_tmp = num_counts[in_row, ]
		denom_tmp = denom_counts[in_row, ]
		
		if (t_paired) {
			num_n = sum(!is.na(num_tmp))
			num_tmp[is.na(num_tmp)] = use_missing
			denom_n = sum(!is.na(denom_tmp))
			denom_tmp[is.na(denom_tmp)] = use_missing
		}	else {
			num_tmp = num_tmp[!is.na(num_tmp)]
			num_n = length(num_tmp)
			if (num_n < 3) {
				num_tmp = rep(use_missing, length(num_samples))
			}
			denom_tmp = denom_tmp[!is.na(denom_tmp)]
			denom_n = length(denom_tmp)
			if (denom_n < 3) {
				denom_tmp = rep(use_missing, length(denom_samples))
			}
			if ((denom_n < 3) & (num_n < 3)) {
				return(NULL)
			}
		}
		t_res = broom::tidy(t.test(num_tmp, denom_tmp, paired = t_paired))
		if (!is.null(paired)) {
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

merge_and_map_metabolomics = function(in_data,
																			chebi_inchikey)
{
	# in_data = tar_read(metabolomics_de_treatment_list)
	# tar_load(chebi_inchikey)
	de_vals = purrr::imap(in_data, \(use_data, id){
		if ("in_ch_i_key" %in% names(use_data)) {
			use_data = use_data |>
				dplyr::mutate(in_chi_key = in_ch_i_key)
		}
		out_de = use_data |>
			dplyr::select(log2FoldChange, p.value, padj, feature_id, in_chi_key) |>
			dplyr::mutate(type = id)
		out_de
	}) |>
		purrr::list_rbind()
	
	de_vals = dplyr::left_join(de_vals, chebi_inchikey, by = c("in_chi_key" = "in_ch_i_key"),
														 relationship = "many-to-many")
	de_vals = de_vals |>
		dplyr::mutate(feature_org = feature_id,
									feature_id = as.character(chebi_id))
	de_vals
}

compare_treatment_patient = function(de_treatment,
																		 de_patient)
{
	# de_treatment = tar_read(metabolomics_de_treatment)
	# de_patient = tar_read(metabolomics_de_patient)
	# 
	# de_treatment = tar_read(rna_de_treatment)
	# de_patient = tar_read(rna_de_patient)
	
	if ("feature_org" %in% names(de_treatment)) {
		de_treatment_mod = de_treatment |>
			dplyr::select(log2FoldChange, p.value, padj, feature_org, type) |>
			dplyr::distinct()
		de_patient_mod = de_patient |>
			dplyr::select(log2FoldChange, p.value, padj, feature_org, type) |>
			dplyr::distinct()
		de_merge = dplyr::left_join(de_treatment_mod, de_patient_mod,
																by = "feature_org",
																suffix = c(".treatment", ".patient"))
	} else {
		de_treatment_mod = de_treatment |>
			dplyr::select(log2FoldChange, pvalue, padj, feature_id) |>
			dplyr::distinct()
		de_patient_mod = de_patient |>
			dplyr::select(log2FoldChange, pvalue, padj, feature_id) |>
			dplyr::distinct()
		de_merge = dplyr::left_join(de_treatment_mod, de_patient_mod,
																by = "feature_id",
																suffix = c(".treatment", ".patient"))
	}
	de_merge
	
}