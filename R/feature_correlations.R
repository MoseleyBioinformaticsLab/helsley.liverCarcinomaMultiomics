feature_correlations = function(rna_collapsed,
																rna_significant,
																bioamines_collapsed,
																lipidomics_collapsed,
																pm_collapsed,
																metabolites_significant,
																matched_samples)
{
	# tar_load(c(rna_collapsed,
	# 				 bioamines_collapsed,
	# 				 lipidomics_collapsed,
	# 				 pm_collapsed,
	# 				 matched_samples))
	# rna_significant = tar_read(rna_de_patient)
	# metabolites_significant = tar_read(metabolomics_de_patient_list)
	rna_sig = rna_significant |>
		dplyr::filter(padj <= 0.05)
	metabolites_sig = metabolites_significant |>
		dplyr::filter(padj <= 0.05)
	
	rna_counts = assays(rna_collapsed)$counts[rna_sig$feature_id, matched_samples]
	bio_counts = assays(bioamines_collapsed)$counts[, matched_samples]
	lipid_counts = assays(lipidomics_collapsed)$counts[, matched_samples]
	pm_counts = assays(pm_collapsed)$counts[, matched_samples]
	
	metabolite_counts = rbind(bio_counts, lipid_counts, pm_counts)
	metabolite_counts = metabolite_counts[metabolites_sig$feature_id, ]
	all_counts = rbind(rna_counts, metabolite_counts)
	
	metabolite_id = metabolites_sig$feature_id
	rna_id = rownames(rna_counts)
	
	use_comparison = tidyr::expand_grid(rna = rna_id, metabolite = metabolite_id)
	
	all_cor = ICIKendallTau::ici_kendalltau(t(all_counts), global_na = 0, include_only = use_comparison, return_matrix = FALSE)
	
	just_cor = all_cor$cor |>
		dplyr::filter(s1 != s2)
	just_cor
}


get_matched_samples = function(rna_collapsed,
															 bioamines_collapsed,
															 lipidomics_collapsed,
															 pm_collapsed)
{
	rna_counts = assays(rna_collapsed)$counts
	bio_counts = assays(bioamines_collapsed)$counts
	lipid_counts = assays(lipidomics_collapsed)$counts
	pm_counts = assays(pm_collapsed)$counts
	
	match_samples = intersect(colnames(rna_counts), 
														intersect(colnames(bio_counts), 
																			intersect(colnames(lipid_counts), colnames(pm_counts))))
	match_samples
}

rna_within_correlation = function(rna_collapsed,
																	rna_significant,
																	matched_samples)
{
	# tar_load(c(rna_collapsed,
	# 		matched_samples))
	# rna_significant = tar_read(rna_de_patient)
	
	rna_sig = rna_significant |>
		dplyr::filter(padj <= 0.05)
	rna_counts = assays(rna_collapsed)$counts[rna_sig$feature_id, matched_samples]
	
	rna_cor = ICIKendallTau::ici_kendalltau(t(rna_counts), global_na = 0, return_matrix = FALSE)
	rna_cor
}

metabolites_within_correlation = function(bioamines_collapsed,
																					lipidomics_collapsed,
																					pm_collapsed,
																					matched_samples)
{
	# tar_load(c(bioamines_collapsed,
	# 							lipidomics_collapsed,
	# 							pm_collapsed,
	# 							matched_samples))
	
	bioamines_counts = assays(bioamines_collapsed)$counts[, matched_samples]
	lipidomics_counts = assays(lipidomics_collapsed)$counts[, matched_samples]
	pm_counts = assays(pm_collapsed)$counts[, matched_samples]
	
	all_counts = rbind(bioamines_counts,
										 lipidomics_counts,
										 pm_counts)
	all_cor = ICIKendallTau::ici_kendalltau(t(all_counts), global_na = 0, return_matrix = FALSE)
	all_cor
}