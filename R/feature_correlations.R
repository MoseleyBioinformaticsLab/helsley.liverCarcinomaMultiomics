feature_correlations = function(rna_collapsed,
																rna_significant,
																bioamines_collapsed,
																lipidomics_collapsed,
																pm_collapsed,
																metabolites_significant,
																matched_samples,
																method = "icikt",
																significant_only = TRUE)
{
	# tar_load(c(rna_collapsed,
	# 				 bioamines_collapsed,
	# 				 lipidomics_collapsed,
	# 				 pm_collapsed,
	# 				 matched_samples))
	# rna_significant = tar_read(rna_de_patient)
	# metabolites_significant = tar_read(metabolomics_de_patient_list)
	# method = "icikt"
	if (significant_only) {
		rna_sig = rna_significant |>
			dplyr::filter(padj <= 0.05)
		metabolites_sig = metabolites_significant |>
			dplyr::filter(padj <= 0.05)
	} else {
		rna_sig = rna_significant
		metabolites_sig = metabolites_significant
	}
	
	
	rna_counts = assays(rna_collapsed)$counts[rna_sig$feature_id, matched_samples]
	bio_counts = assays(bioamines_collapsed)$counts[, matched_samples]
	lipid_counts = assays(lipidomics_collapsed)$counts[, matched_samples]
	pm_counts = assays(pm_collapsed)$counts[, matched_samples]
	
	metabolite_counts = rbind(bio_counts, lipid_counts, pm_counts)
	metabolite_counts = metabolite_counts[metabolites_sig$feature_id, ]
	all_counts = rbind(rna_counts, metabolite_counts)
	
	metabolite_id = metabolites_sig$feature_id
	rna_id = rownames(rna_counts)
	
	# testing this out, comment out later
	# rna_id = rna_id[1:50]
	# metabolite_id = metabolite_id[1:50]
	
	use_comparison = tidyr::expand_grid(rna = rna_id, metabolite = metabolite_id)
	
	out_cor = switch(method,
									 icikt = {
									 	all_cor = ICIKendallTau::ici_kendalltau(t(all_counts), global_na = 0,
									 																					include_only = use_comparison,
									 																					return_matrix = FALSE)
									 	just_cor = all_cor$cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 },
									 spearman = {
									  all_cor = run_correlation(t(all_counts), include_only = use_comparison,
									  													method = "spearman", use = "everything")
									  just_cor = all_cor |>
									  	dplyr::filter(s1 != s2)
									  just_cor
									 },
									 pearson = {
									 	all_counts_na = all_counts
									 	all_counts_na[all_counts == 0] = NA
									 	all_cor = run_correlation(t(all_counts_na), include_only = use_comparison,
									 														method = "pearson", use = "pairwise.complete.obs")
									 	just_cor = all_cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 })
	out_cor$padjust = p.adjust(out_cor$pvalue, method = "BH")
	out_cor$method = method
	out_cor
	return(out_cor)
}

find_rna_metabolite_pairs = function(rna_de_patient,
																		 metabolomics_de_patient_list,
																		 rna_metabolites_all_spearman_sig,
																		 rna_compounds_matrix)
{
	# tar_load(c(rna_de_patient,
	# 					 metabolomics_de_patient_list,
	# 					 rna_metabolites_all_spearman_sig,
	# 					 rna_compounds_matrix))
	
	keep_genes = purrr::map(rna_compounds_matrix[c("lipids", "compounds")], rownames) |> 
		unlist(use.names = FALSE) |>
		unique()
	keep_metabolites = purrr::map(rna_compounds_matrix[c("lipids", "compounds")], colnames) |>
		unlist(use.names = FALSE) |>
		unlist()
	
	trim_spearman_sig = rna_metabolites_all_spearman_sig |>
		dplyr::filter(gene %in% keep_genes, metabolite %in% keep_metabolites)
	
	trim_rna_de_patient = rna_de_patient |>
		dplyr::filter(feature_id %in% keep_genes)
	
	rna_de_spearman = dplyr::inner_join(trim_rna_de_patient, trim_spearman_sig, suffix = c(".de", ".cor"), by = c("feature_id" = "gene"))
	rna_de_spearman = rna_de_spearman |>
		dplyr::mutate(transcript = feature_id)
		
	rna_de_spearman
}

bind_metabolomics_counts = function(bioamines_collapsed,
																			lipidomics_collapsed,
																			pm_collapsed,
																			matched_samples)
{
	# tar_load(c(bioamines_collapsed,
	# 				 lipidomics_collapsed,
	# 				 pm_collapsed,
	# 				 matched_samples))
	
	keep_cols = c("replicate", "sample_id", "treatment", "Treatment")
	keep_rows = c("feature_id", "metabolite_id", "in_chi_key")
	
	bioamines_collapsed = bioamines_collapsed[, matched_samples]
	bioamines_norm = counts(bioamines_collapsed, normalized = TRUE)
	bioamines_rows = rowData(bioamines_collapsed)[, keep_rows]
	bioamines_cols = colData(bioamines_collapsed)[, keep_cols]
	
	lipidomics_collapsed = lipidomics_collapsed[, matched_samples]
	lipidomics_norm = counts(lipidomics_collapsed, normalized = TRUE)
	lipidomics_rows = rowData(lipidomics_collapsed)[, keep_rows]
	lipidomics_cols = colData(lipidomics_collapsed)[, keep_cols]
	
	pm_collapsed = pm_collapsed[, matched_samples]
	pm_norm = counts(pm_collapsed, normalized = TRUE)
	pm_rows = rowData(pm_collapsed)[, keep_rows]
	pm_cols = colData(pm_collapsed)[, keep_cols]
	
	all_norm = rbind(bioamines_norm, lipidomics_norm,
									 pm_norm)
	all_rows = rbind(bioamines_rows, lipidomics_rows, pm_rows)
	all_cols = rbind(bioamines_cols, lipidomics_cols, pm_cols)
	
	all_collapsed = SummarizedExperiment(assays = SimpleList(counts = all_norm), rowData = all_rows, colData = pm_cols)
	
	return(all_collapsed)
}

get_rna_norm_values = function(rna_collapsed,
															 matched_samples)
{
	# tar_load(c(rna_collapsed,
	# 				 matched_samples))
	
	norm_counts = counts(rna_collapsed, normalized = TRUE)
	rna_cols = colData(rna_collapsed)
	rna_ranges = SummarizedExperiment::rowRanges(rna_collapsed)
	rna_norm = SummarizedExperiment(assays = SimpleList(counts = norm_counts), colData = rna_cols,
																	rowRanges = rna_ranges)
	
	return(rna_norm)
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
																	matched_samples,
																	method = "icikt")
{
	# tar_load(c(rna_collapsed,
	# 		matched_samples))
	# rna_significant = tar_read(rna_de_patient)
	
	rna_sig = rna_significant |>
		dplyr::filter(padj <= 0.05)
	rna_counts = assays(rna_collapsed)$counts[rna_sig$feature_id, matched_samples]
	
	out_cor = switch(method,
									 icikt = {
									 	all_cor = ICIKendallTau::ici_kendalltau(t(rna_counts), global_na = 0,
									 																					return_matrix = FALSE)
									 	just_cor = all_cor$cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 },
									 spearman = {
									 	all_cor = run_correlation(t(rna_counts), include_only = NULL,
									 														method = "spearman", use = "everything")
									 	just_cor = all_cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 },
									 pearson = {
									 	all_counts_na = rna_counts
									 	all_counts_na[rna_counts == 0] = NA
									 	all_cor = run_correlation(t(all_counts_na), include_only = NULL,
									 														method = "pearson", use = "pairwise.complete.obs")
									 	just_cor = all_cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 })
	out_cor$padjust = p.adjust(out_cor$pvalue, method = "BH")
	out_cor$method = method
	out_cor
}

metabolites_within_correlation = function(bioamines_collapsed,
																					lipidomics_collapsed,
																					pm_collapsed,
																					matched_samples,
																					method = "icikt")
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
	out_cor = switch(method,
									 icikt = {
									 	all_cor = ICIKendallTau::ici_kendalltau(t(all_counts), global_na = 0,
									 																					return_matrix = FALSE)
									 	just_cor = all_cor$cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 },
									 spearman = {
									 	all_cor = run_correlation(t(all_counts), include_only = NULL,
									 														method = "spearman", use = "everything")
									 	just_cor = all_cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 },
									 pearson = {
									 	all_counts_na = all_counts
									 	all_counts_na[all_counts == 0] = NA
									 	all_cor = run_correlation(t(all_counts_na), include_only = NULL,
									 														method = "pearson", use = "pairwise.complete.obs")
									 	just_cor = all_cor |>
									 		dplyr::filter(s1 != s2)
									 	just_cor
									 })
	out_cor$padjust = p.adjust(out_cor$pvalue, method = "BH")
	out_cor$method = method
	out_cor
	return(out_cor)
}

run_correlation = function(data_matrix, include_only, method, use = "everything")
{
	if (requireNamespace("furrr", quietly = TRUE)) {
		ncore = future::nbrOfWorkers()
		names(ncore) = NULL
		split_fun = furrr::future_map
	} else {
		ncore = 1
		split_fun = purrr::map
	}
	
	n_sample = ncol(data_matrix)
	
	pairwise_comparisons = utils::combn(n_sample, 2)
	named_comparisons = data.frame(s1 = colnames(data_matrix)[pairwise_comparisons[1, ]],
																 s2 = colnames(data_matrix)[pairwise_comparisons[2, ]])
	
	if (!is.null(include_only)) {
		if (is.character(include_only) || is.numeric(include_only)) {
			#message("a vector!")
			# Check each of the comparison vectors against the include_only variable
			# This returns TRUE where they match
			# Use OR to make sure we return everything that should be returned
			s1_include = named_comparisons$s1 %in% include_only
			s2_include = named_comparisons$s2 %in% include_only
			named_comparisons = named_comparisons[(s1_include | s2_include), ]
		} else if (is.list(include_only)) {
			if (length(include_only) == 2) {
				#message("a list!")
				# In this case the include_only is a list of things, so we have to check both
				# of the sets against each of the lists. Again, this returns TRUE where
				# they match. Because we want the things where they
				# are both TRUE (assuming l1[1] goes with l2[1]), we use the AND at the end.
				l1_include = (named_comparisons$s1 %in% include_only[[1]]) | (named_comparisons$s2 %in% include_only[[1]])
				l2_include = (named_comparisons$s1 %in% include_only[[2]]) | (named_comparisons$s2 %in% include_only[[2]])
				named_comparisons = named_comparisons[(l1_include & l2_include), ]
			} else {
				stop("include_only must either be a single vector, or a list of 2 vectors!")
			}
		}
	}
	
	if (nrow(named_comparisons) == 0) {
		stop("nrow(named_comparisons) == 0, did you create include_only correctly?")
	}
	
	n_todo = nrow(named_comparisons)
	n_each = ceiling(n_todo / ncore)
	
	which_core = rep(seq(1, ncore), each = n_each)
	which_core = which_core[1:nrow(named_comparisons)]
	
	named_comparisons$core = which_core
	named_comparisons$cor = Inf
	named_comparisons$pvalue = Inf
	
	
	split_comparisons = split(named_comparisons, named_comparisons$core)
	
	do_split = function(do_comparisons, data_matrix, method, use) {
		#seq_range = seq(in_range[1], in_range[2])
		#print(seq_range)
		
		cor = vector("numeric", nrow(do_comparisons))
		pvalue = cor
		
		for (irow in seq_len(nrow(do_comparisons))) {
			iloc = do_comparisons[irow, 1]
			jloc = do_comparisons[irow, 2]
			
			n_match = sum(!is.na(data_matrix[, iloc]) & !is.na(data_matrix[, jloc]))
			if (n_match >= 4) {
				cor_res = cor.test(data_matrix[, iloc], data_matrix[, jloc], method = method, use = use)
				cor[irow] = cor_res$estimate[[1]]
				pvalue[irow] = cor_res$p.value
			} else {
				cor[irow] = NA
				pvalue[irow] = NA
			}
			
		}
		do_comparisons$cor = cor
		do_comparisons$pvalue = pvalue
		#return(ls())
		do_comparisons
	}
	
	split_cor = split_fun(split_comparisons, do_split, data_matrix, method, use)
	
	all_cor = purrr::list_rbind(split_cor)
	all_cor
}

just_rna_abundances = function(rna_collapsed,
																matched_samples)
{
	rna_counts = assays(rna_collapsed)$counts[, matched_samples]
	rna_counts
}

just_metabolite_abundances = function(bioamines_collapsed,
																			 lipidomics_collapsed,
																			 pm_collapsed,
																			 matched_samples)
{
	bioamines_counts = assays(bioamines_collapsed)$counts[, matched_samples]
	lipidomics_counts = assays(lipidomics_collapsed)$counts[, matched_samples]
	pm_counts = assays(pm_collapsed)$counts[, matched_samples]
	
	all_counts = rbind(bioamines_counts,
										 lipidomics_counts,
										 pm_counts)
	all_counts
}

check_metabolite_correlations = function(metabolites_within_cor,
															metabolomics_de_patient_list)
{
	# tar_load(c(metabolites_within_cor,
	# 					 metabolomics_de_patient_list))
	
	metabolomics_de_patient_list = metabolomics_de_patient_list |>
		dplyr::filter(!is.na(in_chi_key))
	has_multiple = metabolomics_de_patient_list |>
		dplyr::group_by(in_chi_key) |>
		dplyr::summarise(n_metabolite = dplyr::n()) |>
		dplyr::filter(n_metabolite > 1)
	
	multiple_ids = metabolomics_de_patient_list |>
		dplyr::filter(in_chi_key %in% has_multiple$in_chi_key)
	
	split_features = split(multiple_ids$feature_id, multiple_ids$in_chi_key)
	
	correlations = purrr::imap(split_features, \(in_features, id){
		tmp_cor = metabolites_within_cor |>
			dplyr::filter(((s1 %in% in_features[1]) & (s2 %in% in_features[2])) | ((s1 %in% in_features[2]) & (s2 %in% in_features[1])))
		tmp_cor$in_chi_key = id
		tmp_cor
	}) |> 
		purrr::list_rbind()
	correlations
}


find_genes_correlated_lipids = function(metabolomics_enrichment_lipid_binomial,
																				rna_metabolites_all_spearman,
																				metabolomics_de_patient_list,
																				rna_de_patient,
																				binomial_padj = 0.05,
																				cor_padj = 0.01,
																				rna_padj = 0.01)
{
	# tar_load(c(metabolomics_enrichment_lipid_binomial,
	# 					 rna_metabolites_all_spearman,
	# 					 metabolomics_de_patient_list,
	#            rna_de_patient))
	# binomial_padj = 0.05
	# cor_padj = 0.05
	# rna_padj = 0.01
	# 
	# metabolomics_enrichment_lipid_binomial = tar_read(metabolomics_enrichment_reactome_binomial)
	# tar_load(rna_metabolites_all_spearman_sig)
	# binomial_padj = 0.1
	# cor_padj = 0.05
	extra_keep = "class:PS"
	force(binomial_padj)
	force(cor_padj)
	force(rna_padj)
	sig_cor = rna_metabolites_all_spearman |>
		dplyr::filter(padjust <= cor_padj) |>
		dplyr::mutate(transcript = s1, metabolite = s2)
	
	sig_binomial = metabolomics_enrichment_lipid_binomial$stats |>
		dplyr::filter((padjust <= binomial_padj) | (id %in% extra_keep))

	sig_rna = rna_de_patient |>
		dplyr::filter(padj <= rna_padj)
	
	sig_binomial_id = sig_binomial$id
	sig_direction = sig_binomial$direction
	names(sig_direction) = sig_binomial$id
	annotated_lipids = metabolomics_enrichment_lipid_binomial$enrichment@annotation@annotation_features[sig_binomial_id]
	
	metabolomics_de_patient_list = metabolomics_de_patient_list |>
		dplyr::mutate(direction = sign(log2FoldChange))
	
	out_genes = purrr::imap(annotated_lipids, \(lipids, id){
		use_direction = sig_direction[id]
		lipid_direction = metabolomics_de_patient_list |>
			dplyr::filter(feature_id %in% lipids, direction == use_direction) |>
			dplyr::pull(feature_id)
		lipid_cor = sig_cor |>
			dplyr::filter(metabolite %in% lipid_direction, transcript %in% sig_rna$feature_id)
		if (nrow(lipid_cor) > 0) {
			lipid_cor$annotation = id
			return(lipid_cor)
		} else {
			return(NULL)
		}
	})
	list(groups = out_genes,
			 measured = unique(sig_cor$transcript),
			 universe = unique(rna_metabolites_all_spearman$s1))
}

create_rna_compounds_matrix = function(compounds = rna_correlated_interesting_compounds,
																			 lipids = rna_correlated_interesting_lipids,
																			 compound_annotation = metabolomics_feature_list,
																			 rna_annotation = ensembl_uniprot,
																			 all_correlation = rna_metabolites_all_spearman)
{
	# compounds = tar_read(rna_correlated_interesting_compounds)
	# lipids = tar_read(rna_correlated_interesting_lipids)
	# compound_annotation = tar_read(metabolomics_feature_list)
	# rna_annotation = tar_read(ensembl_uniprot)
	# all_correlation = tar_read(rna_metabolites_all_spearman)
	
	all_compounds = purrr::map(compounds$groups, \(in_group){
		in_group |> dplyr::select(transcript, metabolite)
	}) |> 
		purrr::list_rbind()
	
	rna_annotation = rna_annotation |>
		dplyr::transmute(s1 = ensembl_gene_id,
										 symbol = hgnc_symbol) |>
		dplyr::distinct()
	
	compound_annotation = compound_annotation |>
		dplyr::transmute(s2 = feature_id,
										 metabolite = metabolite_id) |>
		dplyr::distinct()
	
	compound_df = all_correlation |>
		dplyr::filter(s1 %in% unique(all_compounds$transcript), s2 %in% unique(all_compounds$metabolite))
	
	compound_matrix = ICIKendallTau::long_df_2_cor_matrix(compound_df, is_square = FALSE)
	
	all_lipids = purrr::map(lipids$groups, \(in_group){
		in_group |> dplyr::select(transcript, metabolite)
	}) |> 
		purrr::list_rbind()
	lipid_df = all_correlation |>
		dplyr::filter(s1 %in% unique(all_lipids$transcript), s2 %in% unique(all_lipids$metabolite))
	
	lipid_matrix = ICIKendallTau::long_df_2_cor_matrix(lipid_df, is_square = FALSE)
	
	compound_labels = compound_annotation |>
		dplyr::transmute(feature_id = s2,
										 label = dplyr::case_when(
										 	nchar(metabolite) > 0 ~ metabolite,
										 	nchar(metabolite) == 0 ~ feature_id
										 ))
	
	rna_labels = rna_annotation |>
		dplyr::transmute(feature_id = s1,
										 label = dplyr::case_when(
										 	nchar(symbol) > 0 ~ symbol,
										 	nchar(symbol) == 0 ~ feature_id
										 ))
	
	all_labels = dplyr::bind_rows(compound_labels, rna_labels)
	
	all_labels = all_labels |>
		dplyr::filter(feature_id %in% c(rownames(lipid_matrix), rownames(compound_matrix), colnames(lipid_matrix), colnames(compound_matrix)))
	
	list(compounds = compound_matrix,
			 lipids = lipid_matrix,
			 labels = all_labels)
	
}

