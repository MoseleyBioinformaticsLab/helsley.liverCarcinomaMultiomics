create_pca_heatmap_figures = function(
	qcqa_figure,
	heatmap_figure,
	color_scales,
	pca_legend_pos = c(0.6, 0.8)
) {
	# qcqa_figure = tar_read(bioamines_qcqa)
	# heatmap_figure = tar_read(bioamines_patient_heatmap)
	# tar_load(color_scales)
	# pca_legend_pos = c(0.6, 0.8)
	use_colors = color_scales$normal_cancer

	hm_grob = draw(heatmap_figure) |> grid.grabExpr()

	pca_figure = refactor_treatment(qcqa_figure) |>
		ggplot(aes(x = PC1, y = PC2, color = Treatment)) +
		geom_polygon(
			stat = "ellipse",
			aes(color = Treatment),
			fill = NA,
			linetype = 2,
			show.legend = FALSE
		) +
		geom_point(size = 2) +
		scale_color_manual(values = use_colors) +
		theme(
			legend.position = "inside",
			legend.position.inside = pca_legend_pos,
			legend.text = element_text(size = rel(0.75))
		) +
		labs(
			x = qcqa_figure$pca_nooutlier_variance$labels[1],
			y = qcqa_figure$pca_nooutlier_variance$labels[2]
		)

	combined_figure = list(pca = pca_figure, heatmap = heatmap_figure)
	return(combined_figure)
}

plot_metabolite_rna = function(
	metabolite_rna_pairs,
	metabolite_collapsed_norm,
	rna_collapsed_norm,
	rna_metabolites_all_spearman_sig,
	color_scales
) {
	# tar_load(c(metabolite_rna_pairs,
	# 				 metabolite_collapsed_norm,
	# 				 rna_collapsed_norm,
	# 				 rna_metabolites_all_spearman_sig,
	# 				 color_scales))

	pair_rows = seq_len(nrow(metabolite_rna_pairs))
	use_samples = base::intersect(
		colnames(metabolite_collapsed_norm),
		colnames(rna_collapsed_norm)
	)
	sample_info = colData(metabolite_collapsed_norm) |> tibble::as_tibble()

	cor_plots = purrr::map(pair_rows, \(in_row) {
		# in_row = 1
		use_transcript = metabolite_rna_pairs$transcript[in_row]
		use_metabolite = metabolite_rna_pairs$metabolite[in_row]

		rna_values = rna_collapsed_norm[use_transcript, use_samples] |>
			dedataset_to_df() |>
			dplyr::mutate(RNA = value)
		metabolite_values = metabolite_collapsed_norm[
			use_metabolite,
			use_samples
		] |>
			dedataset_to_df() |>
			dplyr::mutate(Metabolite = value)

		all_values = dplyr::left_join(
			rna_values,
			metabolite_values[, c("Metabolite", "sample_id", "feature_name")],
			suffix = c(".rna", ".metabolite"),
			by = "sample_id"
		)
		all_values = all_values |>
			dplyr::mutate(
				RNA = log2(RNA + 1),
				Metabolite = log2(Metabolite + 1)
			)
		cor_data = rna_metabolites_all_spearman_sig |>
			dplyr::filter(
				gene %in% use_transcript,
				metabolite %in% use_metabolite
			)

		cor_string = glue::glue(
			"Correlation: {format(cor_data$cor[1], digits = 2)}; Adj-P-Value: {format(cor_data$padjust, digits = 2, scientific = TRUE)}"
		)

		out_plot = all_values |>
			ggplot(aes(x = RNA, y = Metabolite, color = Treatment)) +
			geom_point(size = 2, show.legend = FALSE) +
			scale_color_manual(values = color_scales$normal_cancer) +
			labs(
				x = all_values$feature_name.rna[1],
				y = all_values$feature_name.metabolite[1],
				subtitle = cor_string
			)
		out_plot
	})

	metabolite_rna_pairs$cor_plots = cor_plots
	metabolite_rna_pairs
}

dedataset_to_df = function(dedataset) {
	# dedataset = rna_values
	count_data = assays(dedataset)$counts |> as.matrix()
	meta_data = colData(dedataset) |> tibble::as_tibble()
	out_df = tibble::tibble(
		value = count_data[1, ],
		sample_id = meta_data$sample_id,
		Treatment = meta_data$Treatment
	)
	de_rows = rowData(dedataset) |> tibble::as_tibble()

	if ("name" %in% names(de_rows)) {
		out_df$feature_name = de_rows$name
	} else {
		out_df$feature_name = de_rows$metabolite_id
	}

	out_df
}


create_median_correlation_plots = function(median_cor_list, color_scales) {
	# tar_load(median_cor_list)
	# tar_load(color_scales)

	use_colors = color_scales$normal_cancer

	median_cor_figure = function(in_data, plot_id) {
		# in_data = median_cor_list[[1]]
		# plot_id = names(median_cor_list)[1]

		if (plot_id %in% c("RNA", "Bioamines")) {
			out_sina = in_data |>
				ggplot(aes(
					x = Treatment,
					y = med_cor,
					group = Treatment,
					shape = Outlier,
					color = Treatment
				)) +
				geom_sina(size = 2) +
				scale_color_manual(values = use_colors) +
				guides(color = "none", shape = "none") +
				labs(
					x = "Treatment",
					y = "Median ICI-Kt Correlations",
					subtitle = plot_id
				)
		} else {
			out_sina = in_data |>
				ggplot(aes(
					x = Treatment,
					y = med_cor,
					group = Treatment,
					shape = outlier,
					color = Treatment
				)) +
				geom_sina(size = 2) +
				scale_color_manual(values = use_colors) +
				guides(color = "none") +
				labs(
					x = "Treatment",
					y = "Median ICI-Kt Correlations",
					subtitle = plot_id
				)
		}

		out_sina
	}
	all_sina = purrr::imap(median_cor_list, median_cor_figure)
	all_sina
}

export_median_cor_pptx = function(
	sina_figures = median_correlation_figures,
	ppt_file = "docs/median_correlation_supplemental_figures.pptx"
) {
	# sina_figures = tar_read(median_correlation_figures)
	# ppt_file = "docs/median_correlation_supplemental_figures.pptx"

	new_ppt = officer::read_pptx()
	new_ppt = officer::add_slide(new_ppt, layout = "Blank")

	sina_figures = purrr::map(sina_figures, \(in_figure) {
		in_figure +
			theme(
				axis.title = element_text(size = 9),
				axis.text = element_text(size = 9),
				legend.text = element_text(size = 9),
				legend.title = element_text(size = 10)
			)
	})

	panel_figure = wrap_plots(
		sina_figures,
		ncol = 2,
		nrow = 2,
		tag_level = "keep",
		guides = "collect"
	) +
		plot_annotation(tag_levels = "A")

	panel_dml = rvg::dml(ggobj = panel_figure)
	new_ppt = officer::ph_with(
		new_ppt,
		value = panel_dml,
		location = officer::ph_location(
			left = 1,
			top = 0.25,
			width = 7.5,
			height = 6.5
		)
	)
	figure_caption = officer::fpar(officer::ftext(
		"Sina plots of sample median ICI-Kt correlations to all other samples of the same type. Color indicates disease status, and shape indicates outlier status. A: RNA; B: Bioamines; C: Lipidomics; D: Primary Metabolism.",
		officer::fp_text(color = "black", font.size = 10)
	))
	figure_loc = officer::ph_location(
		left = 1,
		top = 6.5,
		width = 8,
		height = 1
	)
	new_ppt = officer::ph_with(
		new_ppt,
		value = figure_caption,
		location = figure_loc
	)

	print(new_ppt, target = ppt_file)

	NULL
}

export_pca_heatmaps_pptx = function(plot_list, ppt_file) {
	#plot_list = tar_read(pca_heatmap_list)
	#ppt_file = "docs/pca_heatmap_figures.pptx"
	new_ppt = officer::read_pptx()
	for (iplot in plot_list) {
		new_ppt = officer::add_slide(new_ppt, layout = "Blank")
		# grab the PCA part
		pca_dml = rvg::dml(ggobj = iplot$figure$pca)

		# add PCA
		new_ppt = officer::ph_with(
			new_ppt,
			value = pca_dml,
			location = officer::ph_location(
				left = 1,
				top = 1,
				width = 4,
				height = 4
			)
		)

		# write heatmap to tmp file before incorporating it as a PNG
		tmp_file = tempfile("heatmap", fileext = ".png")
		heatmap_tmp = ragg::agg_png(
			tmp_file,
			width = 4,
			height = 4,
			res = 600,
			units = "in",
			background = "white"
		)
		draw(iplot$figure$heatmap)
		dev.off()

		new_ppt = officer::ph_with(
			new_ppt,
			value = officer::external_img(
				tmp_file,
				width = 4,
				height = 4,
				unit = "in"
			),
			location = officer::ph_location(
				left = 5,
				top = 1,
				width = 4,
				height = 4
			)
		)

		figure_caption = officer::fpar(officer::ftext(
			iplot$caption,
			officer::fp_text(color = "black", font.size = 12)
		))
		figure_loc = officer::ph_location(
			left = 1,
			top = 5,
			width = 8,
			height = 1
		)
		new_ppt = officer::ph_with(
			new_ppt,
			value = figure_caption,
			location = figure_loc
		)
	}
	print(new_ppt, target = ppt_file)
	return(ppt_file)
}


refactor_treatment = function(qcqa_obj) {
	# pca_values = rna_qcqa$pca_nooutlier_values
	pca_values = qcqa_obj$pca_nooutlier_values
	pca_values = pca_values |>
		dplyr::mutate(
			Treatment = dplyr::case_when(
				treatment %in% "cancerous" ~ "Cancerous",
				treatment %in% "normal_adjacent" ~ "Normal Adjacent"
			)
		)
	pca_values
}


generate_metabolomics_de_output = function(metabolomics_de_patient_list) {
	data_dictionary = tibble::tribble(
		~header,
		~meaning,
		"baseMean",
		"mean across control samples calculated by DESeq2",
		"log2FoldChange",
		"Fold change of carcinoma / adjacent normal calculated by DESeq2",
		"lfcSE",
		"standard error of fold change",
		"stat",
		"test statistic",
		"pvalue",
		"p-value",
		"padj",
		"adjusted p-value by Benjamini-Hochberg",
		"identifier",
		"WCMC provided identifier",
		"annotation",
		"WCMC annotation",
		"ion species",
		"ionization state of metabolite",
		"in_chi_key",
		"InChIKey for the metabolite",
		"msi_level",
		"measure of how sure WCMC is of the ID",
		"m_z",
		"mass to charge",
		"ret_time_min",
		"retention time from chromatography",
		"esi_mode",
		"whether measured from positive or negative electrospray ionization",
		"feature_id",
		"a proper ID that can be used in R derived from identifier and the type of metabolite",
		"metabolite_id",
		"a column that makes sure if there is a name for this thing, we caught it",
		"type",
		"which method did the metabolite come from",
		"...",
		"mostly columns specific to each metabolite type"
	)
	tab_out = list(
		dictionary = data_dictionary,
		metabolomics = metabolomics_de_patient_list
	)
	tabular_output = openxlsx::write.xlsx(
		tab_out,
		"docs/metabolomics_patient_differential.xlsx",
		overwrite = TRUE
	)
	tabular_output
}

generate_transcriptomics_de_output = function(rna_de_patient) {
	data_dictionary = tibble::tribble(
		~header,
		~meaning,
		"baseMean",
		"mean across control samples calculated by DESeq2",
		"log2FoldChange",
		"Fold change of carcinoma / adjacent normal calculated by DESeq2",
		"lfcSE",
		"standard error of fold change",
		"stat",
		"test statistic",
		"pvalue",
		"p-value",
		"padj",
		"adjusted p-value by Benjamini-Hochberg",
		"feature_id",
		"Ensembl gene ids",
		"name",
		"Gene symbol",
		"biotype",
		"what type of gene is it",
		"description",
		"gene name"
	)
	tab_out = list(dictionary = data_dictionary, metabolomics = rna_de_patient)
	tabular_output = openxlsx::write.xlsx(
		tab_out,
		"docs/transcriptomics_patient_differential.xlsx",
		overwrite = TRUE
	)
	tabular_output
}

generate_abundance_figures = function(
	rna_metabolites_all_spearman_sig,
	rna_abundances,
	metabolite_abundances,
	rna_de_patient,
	color_scale
) {
	tar_load(
		rna_metabolites_all_spearman_sig,
		rna_abundances,
		metabolite_abundances,
		rna_de_patient,
		color_scales
	)
}

generate_correlation_output = function(
	rna_metabolites_all_spearman_sig,
	rna_de_patient,
	metabolomics_de_patient_list,
	lipid_compounds_in_binomial
) {
	# tar_load(c(rna_metabolites_all_spearman_sig,
	# 					 rna_de_patient,
	# 					 metabolomics_de_patient_list,
	# 							lipid_compounds_in_binomial))
	#

	rna_info = rna_de_patient |>
		dplyr::transmute(
			gene = feature_id,
			gene_symbol = name,
			description = description,
			significant_gene = dplyr::case_when(
				padj <= 0.01 ~ "significant",
				padj > 0.01 ~ "not-significant"
			)
		)
	rna_metabolites_all_spearman_sig = dplyr::left_join(
		rna_metabolites_all_spearman_sig,
		rna_info,
		by = "gene",
		relationship = "many-to-many"
	)
	metabolites_info = metabolomics_de_patient_list |>
		dplyr::transmute(
			metabolite = feature_id,
			metabolite_name = metabolite_id,
			metabolite_type = type,
			significant_compound = dplyr::case_when(
				padj <= 0.01 ~ "significant",
				padj > 0.01 ~ "not-significant"
			),
			significant_binomial = dplyr::case_when(
				metabolite %in% lipid_compounds_in_binomial ~ "significant",
				TRUE ~ "not-significant"
			)
		)
	rna_metabolites_all_spearman_sig = dplyr::left_join(
		rna_metabolites_all_spearman_sig,
		metabolites_info,
		by = "metabolite"
	)

	data_dictionary = tibble::tribble(
		~header,
		~meaning,
		"gene",
		"Ensembl gene id",
		"metabolite",
		"internal metabolite feature id",
		"core",
		"what compute core was the correlation calculated with",
		"cor",
		"Spearman correlation value",
		"pvalue",
		"correlation p-value",
		"padjust",
		"correlation adjusted p-value by Benjamini-Hochberg",
		"method",
		"correlation method",
		"gene_symbol",
		"Gene symbol",
		"description",
		"gene name",
		"significant_gene",
		"was the gene considered significant (padjust <= 0.01)",
		"metabolite_name",
		"the metabolite id",
		"metabolite_type",
		"what type of metabolite was it",
		"significant_metabolite",
		"was the metabolite considered significant",
		"significant_binomial",
		"was the metabolite annotated to a significant annotation and in the same direction as a binomial enrichment"
	)
	tab_out = list(
		dictionary = data_dictionary,
		correlation = rna_metabolites_all_spearman_sig
	)
	tabular_output = openxlsx::write.xlsx(
		tab_out,
		"docs/rna_metabolomics_correlations.xlsx",
		overwrite = TRUE
	)
	tabular_output
}

compare_deseq_limma = function(limma_de, deseq_de) {
	limma_de = tar_read(metabolomics_limma_de)
	deseq_de = tar_read(metabolomics_de_patient_list)

	limma_lfc = limma_de |>
		dplyr::transmute(feature_id = feature_id, limma_lfc = logFC)
	deseq_lfc = deseq_de |>
		dplyr::transmute(feature_id = feature_id, deseq2_lfc = log2FoldChange)

	compare_lfc = dplyr::inner_join(limma_lfc, deseq_lfc, by = "feature_id")

	long_lfc = dplyr::bind_rows(
		deseq_lfc |>
			dplyr::transmute(
				feature_id = feature_id,
				lfc = deseq2_lfc,
				source = "deseq2"
			),
		limma_lfc |>
			dplyr::transmute(
				feature_id = feature_id,
				lfc = limma_lfc,
				source = "limma"
			)
	)

	compare_plot = compare_lfc |>
		ggplot(aes(x = limma_lfc, y = deseq2_lfc)) +
		geom_point(alpha = 0.2) +
		geom_abline(slope = 1, color = "red") +
		geom_vline(xintercept = 0, color = "blue") +
		geom_hline(yintercept = 0, color = "blue")

	compare_hist = long_lfc |>
		ggplot(aes(x = lfc, y = source)) +
		geom_vline(xintercept = 0, color = "red") +
		geom_sina(alpha = 0.2)

	ratio_plot = compare_lfc |>
		dplyr::mutate(limma_deseq_ratio = limma_lfc / deseq2_lfc) |>
		dplyr::filter(abs(limma_deseq_ratio) < 20) |>
		ggplot(aes(x = limma_deseq_ratio)) +
		geom_vline(xintercept = 0, color = "red") +
		geom_histogram(bins = 100) +
		scale_y_log10()

	diff_plot = compare_lfc |>
		dplyr::mutate(limma_deseq_diff = limma_lfc - deseq2_lfc) |>
		dplyr::filter(abs(limma_deseq_diff) < 11) |>
		ggplot(aes(x = limma_deseq_diff)) +
		geom_vline(xintercept = 0, color = "red") +
		geom_histogram(bins = 100) +
		scale_y_log10()
}


generate_groups_output = function(
	rna_binomial_interesting_lipids,
	rna_de_patient,
	metabolomics_de_patient_list,
	out_file = "docs/lipid_genes_binomial_groups.xlsx"
) {
	# tar_load(c(rna_binomial_interesting_lipids,
	# 					 rna_de_patient,
	# 					 metabolomics_de_patient_list))
	# out_file = "docs/lipid_genes_binomial_groups.xlsx"

	rna_info = rna_de_patient |>
		dplyr::transmute(
			gene = feature_id,
			gene_symbol = name,
			description = description
		)
	metabolites_info = metabolomics_de_patient_list |>
		dplyr::transmute(
			metabolite = feature_id,
			metabolite_name = metabolite_id,
			metabolite_type = type
		)

	interesting_groups = rna_binomial_interesting_lipids$groups$interesting$groups
	output_groups = purrr::map(interesting_groups, \(in_group) {
		in_group = in_group |>
			dplyr::select(-transcript, -metabolite) |>
			dplyr::rename(gene = s1, metabolite = s2)
		in_group = dplyr::left_join(in_group, rna_info, by = "gene")
		in_group = dplyr::left_join(
			in_group,
			metabolites_info,
			by = "metabolite"
		)
		in_group
	})
	names(output_groups) = gsub(":", "_", names(output_groups))

	data_dictionary = tibble::tribble(
		~header,
		~meaning,
		"gene",
		"Ensembl gene id",
		"metabolite",
		"internal metabolite feature id",
		"core",
		"what compute core was the correlation calculated with",
		"cor",
		"Spearman correlation value",
		"pvalue",
		"correlation p-value",
		"padjust",
		"correlation adjusted p-value by Benjamini-Hochberg",
		"method",
		"correlation method",
		"gene_symbol",
		"Gene symbol",
		"description",
		"gene name",
		"metabolite_name",
		"the metabolite id",
		"metabolite_type",
		"what type of metabolite was it"
	)

	tab_out = c(list(dictionary = data_dictionary), output_groups)
	tabular_output = openxlsx::write.xlsx(tab_out, out_file, overwrite = TRUE)
	tabular_output
}


write_plot_list_includes = function(
	plot_list,
	out_file = "docs/_lipid_binomial_classes.qmd",
	use_id = "binomial_lipid_class_plots"
) {
	load_code = glue::glue(
		"```{{r}}
    plot_list = targets::tar_read({use_id})
    ```"
	)
	plot_code = purrr::map(names(plot_list), \(in_name) {
		id_name = gsub(" ", "-", tolower(in_name))
		glue::glue(
			"```{{r}}
    #| label: fig-{id_name}
    #| fig-cap: Displacement vs {in_name}.
    #| echo: false
    wrap_plots(plot_list[[\"{in_name}\"]], nrow = 1)
    ```"
		)
	}) |>
		purrr::list_c()

	all_code = c(load_code, plot_code)

	cat(
		all_code,
		file = out_file,
		sep = "

",
		append = FALSE
	)
	cli::cli_alert_info(
		"Make sure to have 
  
  {.strong tar_load(include_name)} 
  
  in a code block, and 
  
  {.strong {{{{< include {out_file} >}}}} } 
  
  where you want the figures in the parent file."
	)
	return(list(code = all_code, file = out_file))
}

create_binomial_lipid_overall_class = function(
	binomial_up_down_summary,
	metabolomics_enrichment_lipid_binomial,
	color_scales
) {
	# tar_load(c(binomial_up_down_summary,
	# 					 metabolomics_enrichment_lipid_binomial,
	# 					 color_scales))
	use_colors = color_scales$normal_cancer
	names(use_colors) = c("pos", "neg")

	class_summary = binomial_up_down_summary$class
	binomial_class_labels = metabolomics_enrichment_lipid_binomial$stats |>
		dplyr::filter(grepl("^class\\:", id)) |>
		dplyr::mutate(value = stringr::str_split_i(id, "\\:", 2))
	use_classes = base::intersect(
		class_summary$value,
		binomial_class_labels$value
	)

	class_summary = class_summary |>
		dplyr::filter(value %in% use_classes)
	binomial_class_labels = binomial_class_labels |>
		dplyr::filter(value %in% use_classes) |>
		dplyr::select(value, padjust) |>
		dplyr::mutate(label = format(padjust, digits = 2))

	binomial_class_labels = dplyr::left_join(
		binomial_class_labels,
		class_summary |>
			dplyr::filter(direction_char %in% "pos"),
		by = "value"
	)
	binomial_class_labels = binomial_class_labels |>
		dplyr::mutate(
			label_loc = dplyr::case_when(
				!is.na(n) ~ n + 2,
				TRUE ~ 2
			)
		)

	binomial_class_labels = binomial_class_labels |>
		dplyr::arrange(padjust)

	class_order = binomial_class_labels$value
	binomial_class_labels$value = factor(
		binomial_class_labels$value,
		levels = class_order,
		ordered = TRUE
	)
	class_summary$value = factor(
		class_summary$value,
		levels = class_order,
		ordered = TRUE
	)
	class_summary$direction_char = factor(
		class_summary$direction_char,
		levels = c("pos", "neg"),
		ordered = TRUE
	)

	y_lim = c(-1 * max(abs(class_summary$n)), max(abs(class_summary$n)))
	all_class_plot = class_summary |>
		ggplot(aes(x = value, y = n, fill = direction_char)) +
		scale_fill_manual(values = use_colors) +
		geom_bar(stat = "identity") +
		geom_hline(color = "black", yintercept = 0) +
		geom_text(
			data = binomial_class_labels,
			aes(x = value, y = label_loc, label = label),
			size = 3
		) +
		coord_cartesian(ylim = y_lim) +
		labs(x = "Lipid Class", y = "Decreased in Tumor | Increased in Tumor") +
		theme(legend.position = "none", axis.text.x = element_text(angle = 90))

	return(list(plot = all_class_plot, class_list = class_order))
}


create_binomial_lipid_class_plots = function(
	binomial_up_down_summary,
	metabolomics_enrichment_lipid_binomial,
	color_scales
) {
	# tar_load(c(binomial_up_down_summary,
	# 					 metabolomics_enrichment_lipid_binomial,
	# 					 color_scales))
	use_colors = color_scales$normal_cancer
	names(use_colors) = c("pos", "neg")
	class_summary = binomial_up_down_summary$class
	binomial_class_labels = metabolomics_enrichment_lipid_binomial$stats |>
		dplyr::filter(grepl("^class\\:", id)) |>
		dplyr::mutate(value = stringr::str_split_i(id, "\\:", 2))
	use_classes = base::intersect(
		class_summary$value,
		binomial_class_labels$value
	)

	class_summary = class_summary |>
		dplyr::filter(value %in% use_classes)
	binomial_class_labels = binomial_class_labels |>
		dplyr::filter(value %in% use_classes) |>
		dplyr::select(value, padjust) |>
		dplyr::mutate(label = format(padjust, digits = 2))

	binomial_class_labels = dplyr::left_join(
		binomial_class_labels,
		class_summary |>
			dplyr::filter(direction_char %in% "pos"),
		by = "value"
	)
	binomial_class_labels = binomial_class_labels |>
		dplyr::mutate(
			label_loc = dplyr::case_when(
				!is.na(n) ~ n + 2,
				TRUE ~ 2
			)
		)

	binomial_class_labels = binomial_class_labels |>
		dplyr::arrange(padjust)

	class_order = binomial_class_labels$value
	binomial_class_labels$value = factor(
		binomial_class_labels$value,
		levels = class_order,
		ordered = TRUE
	)
	class_summary$value = factor(
		class_summary$value,
		levels = class_order,
		ordered = TRUE
	)
	class_summary$direction_char = factor(
		class_summary$direction_char,
		levels = c("pos", "neg"),
		ordered = TRUE
	)

	y_lim = c(-1 * max(abs(class_summary$n)), max(abs(class_summary$n)))
	all_class_plot = class_summary |>
		ggplot(aes(x = value, y = n, fill = direction_char)) +
		scale_fill_manual(values = use_colors) +
		geom_bar(stat = "identity") +
		geom_hline(color = "black", yintercept = 0) +
		geom_text(
			data = binomial_class_labels,
			aes(x = value, y = label_loc, label = label),
			size = 3
		) +
		coord_cartesian(ylim = y_lim) +
		labs(x = "Lipid Class", y = "Decreased in Tumor | Increased in Tumor") +
		theme(legend.position = "none", axis.text.x = element_text(angle = 90))

	class_plots = purrr::map(class_order, \(in_class) {
		out_plot = length_db_plot_up_down(
			binomial_up_down_summary$other,
			in_class,
			use_colors
		)
		out_plot[[1]] = out_plot[[1]] + labs(title = in_class)
		out_plot
	})
	names(class_plots) = class_order

	return(class_plots)
}


length_db_plot_up_down = function(
	other_summary,
	which_class = NULL,
	use_colors
) {
	if (!require("ggplot2", quietly = TRUE)) {
		stop("ggplot2 is required to create the plot!")
	}
	if (!is.null(which_class)) {
		other_summary = other_summary |>
			dplyr::filter(class %in% which_class)
	}

	other_summary$type = factor(
		other_summary$type,
		levels = c("total_length", "total_db", "chain_length", "chain_db"),
		ordered = TRUE
	)

	split_type = split(other_summary, other_summary$type)

	out_types = purrr::map(split_type, \(in_type) {
		if (nrow(in_type) > 0) {
			max_val = max(abs(in_type$n))
			in_type$value = as.numeric(in_type$value)
			in_type$direction_char = factor(
				in_type$direction_char,
				levels = c("pos", "neg"),
				ordered = TRUE
			)

			y_lim = c(-1 * max_val, max_val)

			if (grepl("db", in_type$type[1])) {
				x_lab = "# of Double Bonds"
			} else {
				x_lab = "# of Carbons"
			}
			if (grepl("total", in_type$type[1])) {
				sub_title = "Total"
			} else {
				sub_title = "Chain"
			}
			y_lab = NULL

			if (in_type$type[1] %in% "total_length") {
				use_breaks = seq(min(in_type$value), max(in_type$value), 4)
			} else if (in_type$type[1] %in% "total_db") {
				use_breaks = seq(min(in_type$value), max(in_type$value), 4)
			} else if (in_type$type[1] %in% "chain_length") {
				use_breaks = seq(min(in_type$value), max(in_type$value), 2)
			} else if (in_type$type[1] %in% "chain_db") {
				use_breaks = seq(min(in_type$value), max(in_type$value), 2)
			}
			out_plot = in_type %>%
				ggplot(aes(x = value, y = n, fill = direction_char)) +
				scale_fill_manual(values = use_colors) +
				geom_bar(stat = "identity") +
				geom_hline(color = "black", yintercept = 0) +
				coord_cartesian(ylim = y_lim) +
				labs(subtitle = sub_title, x = x_lab, y = y_lab) +
				theme(legend.position = "none")
			return(out_plot)
		} else {
			return(NULL)
		}
	})

	null_plots = purrr::map_lgl(out_types, is.null)
	out_types = out_types[!null_plots]

	return(out_types)
}
