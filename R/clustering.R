label_clusters_first_pass = function(hclust_obj, n_group = 3){
	# hclust_obj = lipid_dist
	# n_group = 3
	out_clusters = cutree(hclust_obj, k = n_group)
	cluster_df = tibble::tibble(feature_id = names(out_clusters),
															cluster = as.character(out_clusters))
	color_labels = ggplot2::scale_color_discrete()$palette(n_group)
	color_df = tibble::tibble(color = color_labels,
														cluster = as.character(seq_len(n_group)))
	cluster_df = dplyr::left_join(cluster_df, color_df, by = "cluster")
	cluster_df
}

label_clusters_second_pass = function(cluster_info, hclust_obj, sub_cluster, n_sub)
{
	# cluster_info = lipid_cluster_first
	# hclust_obj = lipid_dist
	# sub_cluster = 3
	# n_sub = 6
	# 
	sub_cluster = as.character(sub_cluster)
	not_sub = cluster_info |>
		dplyr::filter(!(cluster == sub_cluster)) |>
		dplyr::pull(feature_id)
	pruned_cluster = hclust_obj |>
		dendextend::prune(not_sub)
	sub_color = cluster_info |>
		dplyr::filter(cluster == sub_cluster) |>
		dplyr::pull(color) |>
		unique()
	new_colors = monochromeR::generate_palette(sub_color, modification = "go_both_ways", n_colors = n_sub)
	
	sub_clusters = cutree(pruned_cluster, k = n_sub)
	
	cluster_df = tibble::tibble(feature_id = names(sub_clusters),
															cluster = as.character(sub_clusters))
	color_df = tibble::tibble(color = new_colors,
														cluster = as.character(seq_len(n_sub)))
	cluster_df = dplyr::left_join(cluster_df, color_df, by = "cluster")
	cluster_df = cluster_df |>
		dplyr::mutate(cluster = paste0(sub_cluster, ".", cluster))
	cluster_info$cluster = as.character(cluster_info$cluster)
	
	matches = which(cluster_info$feature_id %in% cluster_df$feature_id)
	
	cluster_info[matches, ] = cluster_df
	
	cluster_info
	
}

create_cluster_color_list = function(cluster_df)
{
	unique_color = cluster_df |>
		dplyr::select(cluster, color) |>
		dplyr::distinct()
	out_colors = unique_color$color
	names(out_colors) = unique_color$cluster
	out_colors
}

#' create class, length, db labels
#' 
#' This function takes the lipid class, total length, and total unsaturations,
#' and creates a label that should be informative for display on a heatmap
create_lipid_labels = function(feature_lipid)
{
	# tar_load(feature_lipid)
	feature_lipid_list = feature_lipid@annotation_features

	lipid_annotation_df = purrr::imap(feature_lipid_list, \(in_feature, in_annotation){
		tibble::tibble(feature_id = in_feature, annotation = in_annotation)
	}) |>
		purrr::list_rbind()

	lipid_class = lipid_annotation_df |>
		dplyr::filter(grepl("^class\\:", annotation)) |>
		dplyr::transmute(feature_id = feature_id,
										 class = gsub("class:", "", annotation))

	lipid_length = lipid_annotation_df |>
		dplyr::filter(grepl("^total_length\\:", annotation)) |>
		dplyr::transmute(feature_id = feature_id,
										 length = gsub("total_length:", "", annotation))
	lipid_db = lipid_annotation_df |>
		dplyr::filter(grepl("^total_db\\:", annotation)) |>
		dplyr::transmute(feature_id = feature_id,
										 db = gsub("total_db:", "", annotation))

	lipid_out_annotation = dplyr::left_join(lipid_class,
											lipid_length, by = "feature_id")
	lipid_out_annotation = dplyr::left_join(lipid_out_annotation,
										lipid_db, by = "feature_id")
		lipid_out_annotation = lipid_out_annotation |>
			dplyr::mutate(annotation_str = paste0(class, "_", length, "_", db))
										
	lipid_out_annotation
}


cluster_create_heatmaps = function(rna_compounds_matrix,
																	 feature_lipid,
																	which = "lipids")
{
	# tar_load(rna_compounds_matrix)	
	# tar_load(feature_lipid)
	# which = "lipids"

	# tar_load(rna_compounds_matrix)	
	# tar_load(feature_lipid)
	# which = "compounds"
	lipid_matrix = rna_compounds_matrix[[which]]
	feature_labels = rna_compounds_matrix$labels
	gene_dist = hclust(dist((1 - lipid_matrix))) |> dendsort::dendsort()
	lipid_dist = hclust(dist(t(1 - lipid_matrix))) |> dendsort::dendsort()

	if (which %in% "lipids") {
		lipid_clusters = label_clusters_first_pass(lipid_dist, 3) |>
			label_clusters_second_pass(lipid_dist, 3, 6) |>
			label_clusters_second_pass(lipid_dist, 2, 2)
	
		gene_clusters = label_clusters_first_pass(gene_dist, 3) |>
			label_clusters_second_pass(gene_dist, 2, 4) |>
			label_clusters_second_pass(gene_dist, 1, 2)
	} else {
		lipid_clusters = label_clusters_first_pass(lipid_dist, 2)
				
		gene_clusters = label_clusters_first_pass(gene_dist, 2) |>
			label_clusters_second_pass(gene_dist, 2, 3) |>
			label_clusters_second_pass(gene_dist, 1, 2)
	}
	
	gene_order = tibble::tibble(feature_id = gene_dist$labels,
															order = order(gene_dist$order))
	lipid_order = tibble::tibble(feature_id = lipid_dist$labels,
																order = order(lipid_dist$order))
	lipid_clusters = dplyr::left_join(lipid_clusters, lipid_order, by = "feature_id")
	gene_clusters = dplyr::left_join(gene_clusters, gene_order, by = "feature_id")
	lipid_colors = list(cluster = create_cluster_color_list(lipid_clusters))
	gene_colors = list(cluster = create_cluster_color_list(gene_clusters))

	full_row_annotation = ComplexHeatmap::HeatmapAnnotation(df = gene_clusters[, "cluster"],
		col = gene_colors, which = "row",
		show_annotation_name = FALSE)
	full_col_annotation = ComplexHeatmap::HeatmapAnnotation(df = lipid_clusters[, "cluster"],
		col = lipid_colors, which = "column")
	cor_map = circlize::colorRamp2(seq(-1, 1, length.out = 20), scico::scico(20, palette = "vik"))
	full_lipid_map = ComplexHeatmap::Heatmap(lipid_matrix, col = cor_map, name = "Spearman",
							top_annotation = full_col_annotation,
							left_annotation = full_row_annotation,
							cluster_rows = gene_dist,
							cluster_columns = lipid_dist,
							show_row_names = FALSE,
							show_column_names = FALSE)
	
	split_gene_clusters = split(gene_clusters, gene_clusters$cluster)

	if (which %in% "lipids") {
		lipid_labels = create_lipid_labels(feature_lipid)
		lipid_chr_labels = lipid_labels$annotation_str
		names(lipid_chr_labels) = lipid_labels$feature_id
	} else {
		lipid_chr_labels = feature_labels$label
		names(lipid_chr_labels) = feature_labels$feature_id
		lipid_chr_labels = lipid_chr_labels[colnames(lipid_matrix)]
	}
	
	gene_lipid_labels = feature_labels$label
	names(gene_lipid_labels) = feature_labels$feature_id
	
	lipid_gene_heatmaps = purrr::imap(split_gene_clusters, \(in_cluster, cluster_id){
		# in_cluster = split_gene_clusters[[1]]
		# cluster_id = names(split_gene_clusters)[1]
		non_cluster = gene_dist$labels[!(gene_dist$labels %in% in_cluster$feature_id)]
		logical_non_cluster = !(gene_dist$labels %in% in_cluster$feature_id)

		tmp_gene_annote = ComplexHeatmap::HeatmapAnnotation(df = in_cluster[, "cluster"],
		col = gene_colors, which = "row",
		show_annotation_name = FALSE)

		subset_dist = dendextend::prune(gene_dist, non_cluster)
		subset_matrix = lipid_matrix[in_cluster$feature_id, ]

		tmp_gene_labels = gene_lipid_labels[rownames(subset_matrix)]
		tmp_lipid_labels = lipid_chr_labels[colnames(subset_matrix)]
		out_heatmap = ComplexHeatmap::Heatmap(subset_matrix, col = cor_map, name = "Cor",
						bottom_annotation = full_col_annotation,
						right_annotation = tmp_gene_annote,
						cluster_rows = subset_dist,
						cluster_columns = lipid_dist,
						show_row_names = TRUE,
						show_column_names = TRUE,
						column_names_side = "bottom",
						column_labels = tmp_lipid_labels,
						row_names_side = "right",
						row_labels = tmp_gene_labels,
						column_names_gp = gpar(fontsize = 6),
						row_names_gp = gpar(fontsize = 6))
		out_heatmap
	})

	out_gene_heatmaps = c(list(full = full_lipid_map),
														lipid_gene_heatmaps)
	if (which %in% "lipids") {
		names(out_gene_heatmaps) = paste0("genes-lipids-", names(out_gene_heatmaps))
	} else {
		names(out_gene_heatmaps) = paste0("genes-compounds-", names(out_gene_heatmaps)) 
	}
	
	return(list(heatmaps = out_gene_heatmaps,
							clusters = list(genes = gene_clusters,
															metabolites = lipid_clusters)))
}

write_lipid_heatmaps_includes = function(plot_list,
	out_file = "docs/_lipid_cluster_maps.qmd",
	use_id = "rna_interesting_lipids_heatmaps")
{
	load_code = glue::glue(
"```{{r}}
#| label: load-lipids
plot_list = targets::tar_read({use_id})
```"
	)
	plot_code = purrr::map(names(plot_list), \(in_name){
		id_name = gsub(" ", "-", tolower(in_name))

		if (grepl("full", in_name)) {
glue::glue("```{{r}}
#| label: fig-lipid-{id_name}
#| fig-cap: Heatmap of lipid cluster {id_name}.
#| fig-height: 20
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		}	else if (grepl("1.1|1.2", in_name)) {
glue::glue("```{{r}}
#| label: fig-lipid-{id_name}
#| fig-cap: Heatmap of lipid cluster {id_name}.
#| fig-height: 15
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		} else {
			glue::glue("```{{r}}
#| label: fig-lipid-{id_name}
#| fig-cap: Heatmap of lipid cluster {id_name}.
#| fig-height: 8
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		}


	}) |> purrr::list_c()

	all_code = c(load_code, plot_code)

	cat(all_code, file = out_file, sep = "\n\n", append = FALSE)
	out_file_nodir = basename(out_file)
	cli::cli_alert_info("Make sure to have 

	{.strong tar_load(include_name)} 

	in a code block, and 

	{.strong {{{{< include {out_file_nodir} >}}}} } 

	where you want the figures in the parent file.")
	return(list(code = all_code, file = out_file))

}

write_compound_heatmaps_includes = function(plot_list,
	out_file = "docs/_compound_cluster_maps.qmd",
	use_id = "rna_interesting_compounds_heatmaps")
{
	load_code = glue::glue(
"```{{r}}
#| label: load-compounds
plot_list = targets::tar_read({use_id})
```"
	)
	plot_code = purrr::map(names(plot_list), \(in_name){
		id_name = gsub(" ", "-", tolower(in_name))

		if (grepl("full", in_name)) {
glue::glue("```{{r}}
#| label: fig-compound-{id_name}
#| fig-cap: Heatmap of compound cluster {id_name}.
#| fig-height: 18
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		}	else if (grepl("1.1", in_name)) {
glue::glue("```{{r}}
#| label: fig-compound-{id_name}
#| fig-cap: Heatmap of compound cluster {id_name}.
#| fig-height: 15
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		} else {
			glue::glue("```{{r}}
#| label: fig-compound-{id_name}
#| fig-cap: Heatmap of compound cluster {id_name}.
#| fig-height: 8
#| fig-width: 15
#| echo: false
plot_list[[\"{in_name}\"]]
```")
		}


	}) |> purrr::list_c()

	all_code = c(load_code, plot_code)

	cat(all_code, file = out_file, sep = "\n\n", append = FALSE)
	out_file_nodir = basename(out_file)
	cli::cli_alert_info("Make sure to have 

	{.strong tar_load(include_name)} 

	in a code block, and 

	{.strong {{{{< include {out_file_nodir} >}}}} } 

	where you want the figures in the parent file.")
	return(list(code = all_code, file = out_file))

}

cluster_add_feature_info = function(clusters,
																		rna_list,
																		metabolites_list,
																		feature_lipid)
{
	# clusters = tar_read(rna_interesting_lipids_hc)$clusters
	# rna_list = tar_read(rna_de_patient)
	# metabolites_list = tar_read(metabolomics_de_patient_list)
	# tar_load(feature_lipid)

	gene_data = dplyr::left_join(clusters$genes, rna_list, by = "feature_id") |>
		dplyr::arrange(order)
	metabolites_data = dplyr::left_join(clusters$metabolites, metabolites_list, by = "feature_id") |>
		dplyr::arrange(order)

	frac_lipids = sum(grepl("lipids", metabolites_data$feature_id)) / nrow(metabolites_data)

	if (frac_lipids > 0.5) {
		lipid_info = create_lipid_labels(feature_lipid)
		metabolites_data = dplyr::left_join(metabolites_data, lipid_info, by = "feature_id")
	}

	return(list(genes = gene_data,
				 metabolites = metabolites_data))
}