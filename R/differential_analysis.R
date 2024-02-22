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
	# rna_se = tar_read(rna_paired)
	# which = "patient"
	if (which %in% "treatment") {
		design(rna_se) = ~ treatment
	} else {
		design(rna_se) = ~ patient + treatment
		
	}
	
	rna_deseq = DESeq(rna_se)
	rna_results = results(rna_deseq, contrast = contrast)
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
	contrast_var = contrast[1]
	numerator = contrast[2]
	denominator = contrast[3]
	
	num_samples = data_se[[contrast_var]] == numerator
	denom_samples = data_se[[contrast_var]] == denominator
	
	if (!is.null(paired)) {
		# do something with the pairing variable
	} else {
		counts = assays(data_se)$counts
		counts[]
	}
	NULL
}