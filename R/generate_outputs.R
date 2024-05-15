generate_metabolomics_de_output = function(metabolomics_de_patient_list)
{
	data_dictionary = tibble::tribble(
		~header, ~meaning,
		"baseMean", "mean across control samples calculated by DESeq2",
		"log2FoldChange", "Fold change of carcinoma / adjacent normal calculated by DESeq2",
		"lfcSE", "standard error of fold change",
		"stat", "test statistic",
		"pvalue", "p-value",
		"padj", "adjusted p-value by Benjamini-Hochberg",
		"identifier", "WCMC provided identifier",
		"annotation", "WCMC annotation",
		"ion species", "ionization state of metabolite",
		"in_chi_key", "InChIKey for the metabolite",
		"msi_level", "measure of how sure WCMC is of the ID",
		"m_z", "mass to charge",
		"ret_time_min", "retention time from chromatography",
		"esi_mode", "whether measured from positive or negative electrospray ionization",
		"feature_id", "a proper ID that can be used in R derived from identifier and the type of metabolite",
		"metabolite_id", "a column that makes sure if there is a name for this thing, we caught it",
		"type", "which method did the metabolite come from",
		"...", "mostly columns specific to each metabolite type"
	)
	tab_out = list(dictionary = data_dictionary,
								 metabolomics = metabolomics_de_patient_list)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				"docs/metabolomics_patient_differential.xlsx",
																				overwrite = TRUE)
	tabular_output
}

generate_transcriptomics_de_output = function(rna_de_patient)
{
	data_dictionary = tibble::tribble(
		~header, ~meaning,
		"baseMean", "mean across control samples calculated by DESeq2",
		"log2FoldChange", "Fold change of carcinoma / adjacent normal calculated by DESeq2",
		"lfcSE", "standard error of fold change",
		"stat", "test statistic",
		"pvalue", "p-value",
		"padj", "adjusted p-value by Benjamini-Hochberg",
		"feature_id", "Ensembl gene ids",
		"name", "Gene symbol",
		"biotype", "what type of gene is it",
		"description", "gene name"
	)
	tab_out = list(dictionary = data_dictionary,
								 metabolomics = rna_de_patient)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				"docs/transcriptomics_patient_differential.xlsx",
																				overwrite = TRUE)
	tabular_output
}

generate_correlation_output = function(rna_metabolites_all_spearman,
																			rna_de_patient,
																			metabolomics_de_patient_list)
{
	tar_load(c(rna_metabolites_all_spearman,
						 rna_de_patient,
						 metabolomics_de_patient_list))
	
	
	rna_metabolites_all_spearman = rna_metabolites_all_spearman |>
		dplyr::rename(gene = s1, metabolite = s2)
	rna_info = rna_de_patient |>
		dplyr::transmute(gene = feature_id, gene_symbol = name, description = description)
	rna_metabolites_all_spearman = dplyr::left_join(rna_metabolites_all_spearman, rna_info, by = "gene",
																									relationship = "many-to-many")
	metabolites_info = metabolomics_de_patient_list |>
		dplyr::transmute(metabolite = feature_id, metabolite_name = metabolite_id, metabolite_type = type)
	rna_metabolites_all_spearman = dplyr::left_join(rna_metabolites_all_spearman, metabolites_info, by = "metabolite")
	
	data_dictionary = tibble::tribble(
		~header, ~meaning,
		"gene", "Ensembl gene id",
		"metabolite", "internal metabolite feature id",
		"core", "what compute core was the correlation calculated with",
		"cor", "Spearman correlation value",
		"pvalue", "correlation p-value",
		"padjust", "correlation adjusted p-value by Benjamini-Hochberg",
		"method", "correlation method",
		"gene_symbol", "Gene symbol",
		"description", "gene name",
		"metabolite_name", "the metabolite id",
		"metabolite_type", "what type of metabolite was it"
	)
	tab_out = list(dictionary = data_dictionary,
								 correlation = rna_metabolites_all_spearman)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				"docs/rna_metabolomics_correlations.xlsx",
																				overwrite = TRUE)
	tabular_output
}

generate_groups_output = function()
{
	
}