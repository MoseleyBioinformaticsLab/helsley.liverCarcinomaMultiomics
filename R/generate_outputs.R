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

generate_correlation_output = function(rna_metabolites_all_spearman_sig,
																			rna_de_patient,
																			metabolomics_de_patient_list)
{
	# tar_load(c(rna_metabolites_all_spearman,
	# 					 rna_de_patient,
	# 					 metabolomics_de_patient_list))
	# 
	
	
	rna_info = rna_de_patient |>
		dplyr::transmute(gene = feature_id, gene_symbol = name, description = description)
	rna_metabolites_all_spearman_sig = dplyr::left_join(rna_metabolites_all_spearman_sig, rna_info, by = "gene",
																									relationship = "many-to-many")
	metabolites_info = metabolomics_de_patient_list |>
		dplyr::transmute(metabolite = feature_id, metabolite_name = metabolite_id, metabolite_type = type)
	rna_metabolites_all_spearman_sig = dplyr::left_join(rna_metabolites_all_spearman_sig, metabolites_info, by = "metabolite")
	
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
								 correlation = rna_metabolites_all_spearman_sig)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				"docs/rna_metabolomics_correlations.xlsx",
																				overwrite = TRUE)
	tabular_output
}

generate_groups_output = function(rna_binomial_interesting_lipids,
																	rna_de_patient,
																	metabolomics_de_patient_list,
																	out_file = "docs/lipid_genes_binomial_groups.xlsx")
{
	# tar_load(c(rna_binomial_interesting_lipids,
	# 					 rna_de_patient,
	# 					 metabolomics_de_patient_list))
	# out_file = "docs/lipid_genes_binomial_groups.xlsx"
	
	rna_info = rna_de_patient |>
		dplyr::transmute(gene = feature_id, gene_symbol = name, description = description)
	metabolites_info = metabolomics_de_patient_list |>
		dplyr::transmute(metabolite = feature_id, metabolite_name = metabolite_id, metabolite_type = type)
	
	interesting_groups = rna_binomial_interesting_lipids$groups$interesting$groups
	output_groups = purrr::map(interesting_groups, \(in_group){
		in_group = in_group |>
			dplyr::select(-transcript, -metabolite) |>
			dplyr::rename(gene = s1, metabolite = s2) 
		in_group = dplyr::left_join(in_group, rna_info, by = "gene")
		in_group = dplyr::left_join(in_group, metabolites_info, by = "metabolite")
		in_group
	})
	names(output_groups) = gsub(":", "_", names(output_groups))
	
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
	
	tab_out = c(list(dictionary = data_dictionary),
							output_groups)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				out_file,
																				overwrite = TRUE)
	tabular_output
}


write_plot_list_includes = function(plot_list,
																		out_file = "docs/_lipid_binomial_classes.qmd",
																		use_id = "binomial_lipid_class_plots")
{
	load_code = glue::glue(
		"```{{r}}
    plot_list = targets::tar_read({use_id})
    ```"
	)
	plot_code = purrr::map(names(plot_list), \(in_name){
		id_name = gsub(" ", "-", tolower(in_name))
		glue::glue("```{{r}}
    #| label: fig-{id_name}
    #| fig-cap: Displacement vs {in_name}.
    #| echo: false
    wrap_plots(plot_list[[\"{in_name}\"]], nrow = 1)
    ```")
		
	}) |> purrr::list_c()
	
	all_code = c(load_code, plot_code)
	
	cat(all_code, file = out_file, sep = "

", append = FALSE)
	cli::cli_alert_info("Make sure to have 
  
  {.strong tar_load(include_name)} 
  
  in a code block, and 
  
  {.strong {{{{< include {out_file} >}}}} } 
  
  where you want the figures in the parent file.")
	return(list(code = all_code, file = out_file))
	
}

create_binomial_lipid_overall_class = function(binomial_up_down_summary,
																	metabolomics_enrichment_lipid_binomial)
{
	# tar_load(c(binomial_up_down_summary,
	# 					 metabolomics_enrichment_lipid_binomial))
	
	class_summary = binomial_up_down_summary$class
	binomial_class_labels = metabolomics_enrichment_lipid_binomial$stats |>
		dplyr::filter(grepl("^class\\:", id)) |>
		dplyr::mutate(value = stringr::str_split_i(id, "\\:", 2))
	use_classes = base::intersect(class_summary$value, binomial_class_labels$value)
	
	class_summary = class_summary |>
		dplyr::filter(value %in% use_classes)
	binomial_class_labels = binomial_class_labels |>
		dplyr::filter(value %in% use_classes) |>
		dplyr::select(value, padjust) |>
		dplyr::mutate(label = format(padjust, digits = 2))
	
	binomial_class_labels = dplyr::left_join(binomial_class_labels, class_summary |>
																					 	dplyr::filter(direction_char %in% "pos"),
																					 by = "value")
	binomial_class_labels = binomial_class_labels |>
		dplyr::mutate(label_loc = dplyr::case_when(
			!is.na(n) ~ n + 2,
			TRUE ~ 2))
	
	binomial_class_labels = binomial_class_labels |>
		dplyr::arrange(padjust)
	
	class_order = binomial_class_labels$value
	binomial_class_labels$value = factor(binomial_class_labels$value, levels = class_order, ordered = TRUE)
	class_summary$value = factor(class_summary$value, levels = class_order, ordered = TRUE)
	class_summary$direction_char = factor(class_summary$direction_char, levels = c("pos", "neg"), ordered = TRUE)
	
	y_lim = c(-1 * max(abs(class_summary$n)), max(abs(class_summary$n)))
	all_class_plot = class_summary |>
		ggplot(aes(x = value, y = n, fill = direction_char)) +
		scale_fill_discrete() +
		geom_bar(stat = "identity") +
		geom_hline(color = "black", yintercept = 0) +
		geom_text(data = binomial_class_labels, aes(x = value, y = label_loc, label = label), size = 3) +
		coord_cartesian(ylim = y_lim) +
		labs(x = "Lipid Class", y = "Downchanged / Upchanged") +
		theme(legend.position = "none", axis.text.x = element_text(angle = 90))
	
	return(list(plot = all_class_plot,
			 class_list = class_order))
}


create_binomial_lipid_class_plots = function(binomial_up_down_summary,
																						 metabolomics_enrichment_lipid_binomial)
{
	# tar_load(c(binomial_up_down_summary,
	# 					 metabolomics_enrichment_lipid_binomial))
	
	class_summary = binomial_up_down_summary$class
	binomial_class_labels = metabolomics_enrichment_lipid_binomial$stats |>
		dplyr::filter(grepl("^class\\:", id)) |>
		dplyr::mutate(value = stringr::str_split_i(id, "\\:", 2))
	use_classes = base::intersect(class_summary$value, binomial_class_labels$value)
	
	class_summary = class_summary |>
		dplyr::filter(value %in% use_classes)
	binomial_class_labels = binomial_class_labels |>
		dplyr::filter(value %in% use_classes) |>
		dplyr::select(value, padjust) |>
		dplyr::mutate(label = format(padjust, digits = 2))
	
	binomial_class_labels = dplyr::left_join(binomial_class_labels, class_summary |>
																					 	dplyr::filter(direction_char %in% "pos"),
																					 by = "value")
	binomial_class_labels = binomial_class_labels |>
		dplyr::mutate(label_loc = dplyr::case_when(
			!is.na(n) ~ n + 2,
			TRUE ~ 2))
	
	binomial_class_labels = binomial_class_labels |>
		dplyr::arrange(padjust)
	
	class_order = binomial_class_labels$value
	binomial_class_labels$value = factor(binomial_class_labels$value, levels = class_order, ordered = TRUE)
	class_summary$value = factor(class_summary$value, levels = class_order, ordered = TRUE)
	class_summary$direction_char = factor(class_summary$direction_char, levels = c("pos", "neg"), ordered = TRUE)
	
	y_lim = c(-1 * max(abs(class_summary$n)), max(abs(class_summary$n)))
	all_class_plot = class_summary |>
		ggplot(aes(x = value, y = n, fill = direction_char)) +
		scale_fill_discrete() +
		geom_bar(stat = "identity") +
		geom_hline(color = "black", yintercept = 0) +
		geom_text(data = binomial_class_labels, aes(x = value, y = label_loc, label = label), size = 3) +
		coord_cartesian(ylim = y_lim) +
		labs(x = "Lipid Class", y = "Downchanged / Upchanged") +
		theme(legend.position = "none", axis.text.x = element_text(angle = 90))
	
	class_plots = purrr::map(class_order, \(in_class){
		out_plot = mu_plot_up_down_length_db(binomial_up_down_summary$other, in_class)
		out_plot[[1]] = out_plot[[1]] + labs(title = in_class)
		out_plot
	})
	names(class_plots) = class_order
	
	return(class_plots)
}