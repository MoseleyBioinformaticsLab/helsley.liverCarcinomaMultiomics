get_ensembl_uniprot = function(version = "111")
{
	# version = "111"
	ensembl = useEnsembl(biomart = "ensembl", 
											 dataset = "hsapiens_gene_ensembl", 
											 version = version)
	
	use_attributes = c("ensembl_gene_id",
										 "hgnc_symbol",
										 "uniprot_gn_symbol",
										 "uniprot_gn_id",
										 "description")
	feature_data = getBM(attributes = use_attributes,
											 mart = ensembl)
	
	feature_data
}

get_ensembl_entrez = function(version = "111")
{
	# version = "111"
	ensembl = useEnsembl(biomart = "ensembl", 
											 dataset = "hsapiens_gene_ensembl", 
											 version = version)
	
	use_attributes = c("ensembl_gene_id",
										 "hgnc_symbol",
										 "entrezgene_id",
										 "description")
	feature_data = getBM(attributes = use_attributes,
											 mart = ensembl)
	
	feature_data
}


create_reactome_gene_annotations = function(reactome_file,
																			 target_species = "Homo sapiens")
{
	# tar_load(reactome_file)
	# target_species = "Homo sapiens"
	all_reactome = readr::read_tsv(reactome_file, col_names = FALSE)
	names(all_reactome) = c("ensembl", "pathway", "link", "description", "evidence", "species")
	all_reactome = all_reactome |>
		dplyr::filter(species %in% target_species) |>
		dplyr::filter(grepl("^ENSG", ensembl))
	
	split_reactome = split(all_reactome$ensembl, all_reactome$pathway)
	split_reactome = purrr::map(split_reactome, unique)
	
	reactome_descriptions = all_reactome$description
	names(reactome_descriptions) = all_reactome$pathway
	reactome_descriptions = reactome_descriptions[names(split_reactome)]
	
	out_annotation = categoryCompare2::annotation(annotation_features = split_reactome,
																								annotation_type = "reactome",
																								description = reactome_descriptions,
																								feature_type = "ensembl_gene")	
	return(out_annotation)
}

create_reactome_chebi_annotations = function(chebi_reactome_file,
																						 target_species = "Homo sapiens")
{
	# tar_load(chebi_reactome_file)
	# target_species = "Homo sapiens"
	all_reactome = readr::read_tsv(chebi_reactome_file, col_names = FALSE)
	names(all_reactome) = c("chebi", "pathway", "link", "description", "evidence", "species")
	all_reactome$chebi = as.character(all_reactome$chebi)
	all_reactome = all_reactome |>
		dplyr::filter(species %in% target_species)
	
	split_reactome = split(all_reactome$chebi, all_reactome$pathway)
	split_reactome = purrr::map(split_reactome, unique)
	
	reactome_descriptions = all_reactome$description
	names(reactome_descriptions) = all_reactome$pathway
	reactome_descriptions = reactome_descriptions[names(split_reactome)]
	
	out_annotation = categoryCompare2::annotation(annotation_features = split_reactome,
																								annotation_type = "reactome",
																								description = reactome_descriptions,
																								feature_type = "chebi")	
	return(out_annotation)
}


create_go_annotations = function(ensembl_uniprot,
																 go_file,
																 namespace_file)
{
	# tar_load(ensembl_uniprot)
	# tar_load(go_file)
	# tar_load(namespace_file)
	feature_translation = ensembl_uniprot |>
		dplyr::transmute(from = uniprot_gn_id,
										 to = ensembl_gene_id) |>
		dplyr::filter(nchar(from) > 0, nchar(to) > 0) |>
		dplyr::distinct()
	
	out_annotation = gocats_to_annotation(go_file,
																				namespace_file,
																				feature_type = "ensembl_gene",
																				feature_translation = feature_translation)
	out_annotation
}

group_annotations = function(group_enrichment,
															similarity_cutoff = 0.8)
{
	# group_enrichment = tar_read(rna_treatment_enrichment_pn_go)
	# similarity_cutoff = 0.8
	just_enrich = group_enrichment$enrich
	
	enrich_graph = generate_annotation_graph(just_enrich)
	enrich_graph = suppressMessages(remove_edges(enrich_graph, similarity_cutoff))
	enrich_assign = annotation_combinations(enrich_graph)
	enrich_assign = assign_colors(enrich_assign)
	
	enrich_communities = assign_communities(enrich_graph)
	enrich_comm_labels = label_communities(enrich_communities, just_enrich@annotation)
	
	enrich_table = table_from_graph(enrich_graph, enrich_assign, enrich_comm_labels)
	enrich_table
}

group_annotations_each = function(group_enrichment,
														 similarity_cutoff = 0.8)
{
	# group_enrichment = tar_read(rna_treatment_enrichment_go)
	# similarity_cutoff = 0.8
	just_enrich = group_enrichment$enrich
	just_stats = group_enrichment$stats
	
	enrich_graph = generate_annotation_graph(just_enrich)
	enrich_graph = remove_edges(enrich_graph, similarity_cutoff)
	enrich_assign = annotation_combinations(enrich_graph)
	enrich_assign = assign_colors(enrich_assign)
	
	enrich_communities = assign_communities(enrich_graph)
	enrich_comm_labels = label_communities(enrich_communities, just_enrich@annotation)
	
	sig_each = just_enrich@statistics@significant@significant
	
	each_table = purrr::imap(just_stats, \(in_stat, stat_id){
		# in_stat = just_stats[[1]]
		# stat_id = names(just_stats)[1]
		sig_tmp = rownames(sig_each)[sig_each[, stat_id]]
		tmp_table = in_stat[in_stat$ID %in% sig_tmp, ]
		
		tmp_table
	})
	
	enrich_comm_df = purrr::map(enrich_comm_labels, \(in_label){
		data.frame(group = in_label$label, ID = in_label$members)
	}) |>
		purrr::list_rbind()
	
	each_table_group = purrr::map(each_table, \(in_table){
		out_table = dplyr::left_join(in_table, enrich_comm_df, by = "ID")
		out_table = out_table |>
			dplyr::mutate(group = dplyr::case_when(
				is.na(group) ~ "",
				TRUE ~ group
			))
		out_table
	})
	each_table_group
}


run_enrichment = function(de_values,
													annotation_obj,
													padj_cutoff = 0.01,
													keep_group = NULL)
{
	# de_values = tar_read(rna_de_treatment)
	# annotation_obj = tar_read(ensembl_go)
	# padj_cutoff = 0.01
	de_entries = de_values |>
		dplyr::filter(padj <= padj_cutoff)
	
	de_all = de_entries$feature_id
	de_neg = de_entries |>
		dplyr::filter(log2FoldChange <= 0) |>
		dplyr::pull(feature_id)
	de_pos = de_entries |>
		dplyr::filter(log2FoldChange > 0) |>
		dplyr::pull(feature_id)
	
	universe = de_values$feature_id
	
	enrich_all = hypergeometric_feature_enrichment(
		new("hypergeom_features", significant = de_all,
				universe = universe, annotation = annotation_obj),
		p_adjust = "BH"
	)
	
	stats_all = extract_enrich_stats(enrich_all)
	enrich_neg = hypergeometric_feature_enrichment(
		new("hypergeom_features", significant = de_neg,
				universe = universe, annotation = annotation_obj),
		p_adjust = "BH"
	)
	stats_neg = extract_enrich_stats(enrich_neg)
	enrich_pos = hypergeometric_feature_enrichment(
		new("hypergeom_features", significant = de_pos,
				universe = universe, annotation = annotation_obj),
		p_adjust = "BH"
	)
	stats_pos = extract_enrich_stats(enrich_pos)
	enrich_list = list(all = enrich_all,
										 pos = enrich_pos,
										 neg = enrich_neg)
	stats_list = list(all = stats_all,
										pos = stats_pos,
										neg = stats_neg)
	if (!is.null(keep_group)) {
		enrich_list = enrich_list[keep_group]
		stats_list = stats_list[keep_group]
	}
	enrich_comb = combine_enrichments(enrich_list)
	enrich_sig = get_significant_annotations(enrich_comb, padjust <= 0.01, counts >= 2, counts <= 500)
	
	return(list(enrich = enrich_sig,
							stats = stats_list))
}

map_chebi_to_inchikey = function(inchi_df,
																 inchi_dir = "data/chebi",
																 inchikey_hash)
{
	# tar_load(inchi_df)
	# inchi_dir = "data/chebi"
	# tar_load(inchikey_hash)

	key_files = file.path(inchi_dir, dir(inchi_dir, pattern = ".inchikey"))
	
	key_values = purrr::map(key_files, function(in_file){
		chebi_id = gsub(".inchikey", "", basename(in_file))
		key = readLines(in_file)
		if (length(key) == 0) {
			return(NULL)
		}
		data.frame(chebi_id = chebi_id, in_ch_i_key = key)
	}) |>
		purrr::list_rbind()
	key_values$chebi_id = as.numeric(key_values$chebi_id)
	
	key_map = dplyr::left_join(key_values, inchi_df, by = "chebi_id")
	key_map
}

write_inchi = function(chebi_inchi_file,
											 inchi_dir = "data/chebi")
{
	# tar_load(chebi_inchi_file)

	chebi_inchi = readr::read_table(chebi_inchi_file) |>
		janitor::clean_names() |>
		dplyr::distinct()
	if (sum(duplicated(chebi_inchi$chebi_id)) > 0) {
		stop("There are duplicate ChEBI IDs!")
	}
	chebi_list = split(chebi_inchi$in_ch_i, chebi_inchi$chebi_id)
	purrr::iwalk(chebi_list, \(in_chi, chebi){
		out_file = file.path(inchi_dir, paste0(chebi, ".inchi"))
		cat(in_chi, file = out_file, sep = "\n")
	})
	chebi_inchi
}

get_inchikey_hash = function(inchi_df,
														 inchikey_directory = "data/chebi")
{
	all_files = paste(sort(dir(inchikey_directory, pattern = "*.inchi*")), collapse = ":")
	if (!grepl("inchikey", all_files) || (nchar(all_files) == 0)) {
		stop("No inchikey files exist. Need to run Open Babel first!")
	}
	digest::digest(all_files, algo = "sha256")
}

# to convert the inchi files to inchikey files
# 
# obabel [1]*.inchi -oinchikey -m -e
# obabel [2]*.inchi -oinchikey -m -e
# obabel [3]*.inchi -oinchikey -m -e
# ...
# ...
# obabel [9]*.inchi -oinchikey -m -e
# 

write_goeach_to_excel = function(go_stuff,
																 reactome_stuff)
{
	# go_stuff = tar_read(rna_patient_enrichment_grouped_eachgo)
	# reactome_stuff = tar_read(rna_patient_enrichment_grouped_eachreactome)
	data_dictionary = tibble::tribble(
		~table, ~header, ~meaning,
		"GO*", "p", "hypergeometric p-value",
		"GO*", "odds", "hypergeometric odds ratio, roughly counts / expected",
		"GO*", "expected", "expected number of genes based on number differential",
		"GO*", "counts", "number of genes observed with that annotation in differential",
		"GO*", "padjust", "Benjamini-Hochberg adjusted p-value",
		"GO*", "ID", "identifier of the annotation",
		"GO*", "description", "text description of the annotation",
		"GO*", "group", "is the annotation part of a group of highly related annotations",
		"Reactome*", "p", "hypergeometric p-value",
		"Reactome*", "odds", "hypergeometric odds ratio, roughly counts / expected",
		"Reactome*", "expected", "expected number of genes based on number differential",
		"Reactome*", "counts", "number of genes observed with that annotation in differential",
		"Reactome*", "padjust", "Benjamini-Hochberg adjusted p-value",
		"Reactome*", "ID", "identifier of the annotation",
		"Reactome*", "description", "text description of the annotation",
		"Reactome*", "group", "is the annotation part of a group of highly related annotations")
	
	names(go_stuff) = paste0("GO-", names(go_stuff))
	names(reactome_stuff) = paste0("Reactome-", names(reactome_stuff))
	all_out = c(go_stuff, reactome_stuff)
	all_out = purrr::map(all_out, \(x){
		x = x |>
			dplyr::arrange(dplyr::desc(group), padjust)
		x
	})
	tab_out = c(list(dictionary = data_dictionary),
							all_out)
	tabular_output = openxlsx::write.xlsx(tab_out,
																				"docs/transcriptomics_enrichments.xlsx",
																				overwrite = TRUE)
	tabular_output
}

find_interesting_gm_groups = function(rna_metabolites_spearman,
																			rna_within_cor,
																			metabolites_within_cor,
																			rna_patient_enrichment_grouped_eachgo,
																			rna_patient_enrichment_go,
																			cor_cutoff = 0.01,
																			direction = "neg")
{
	# tar_load(c(
	# 	rna_metabolites_spearman,
	# 	rna_within_cor,
	# 	metabolites_within_cor,
	# 	rna_patient_enrichment_grouped_eachgo,
	# 	rna_patient_enrichment_go
	# ))
	# cor_cutoff = 0.01
	# direction = "pos"
	
	
	
	force(cor_cutoff)
	match_str = "MF:lipoprotein particle binding|CC:vesicle lumen|BP:regulation of tube size|BP:regulation of lipase activity|lipo.*|CC:protein-lipid complex|BP:regulation of spindle checkpoint|BP:regulation of chromosome segregation|BP:negative regulation of cell cycle"
	
	sig_cor = rna_metabolites_spearman |>
		dplyr::filter(padjust <= cor_cutoff) |>
		dplyr::transmute(gene = s1, metabolite = s2, cor = cor, pvalue = pvalue,
										 padjust = padjust)
	cor_genes = unique(sig_cor$gene)
	
	use_go = rna_patient_enrichment_grouped_eachgo[[direction]]
	interesting_go = use_go |>
		dplyr::filter(grepl(match_str, group))
	
	split_go = split(interesting_go$ID, interesting_go$group)
	
	rna_sig = rna_patient_enrichment_go$enrich@enriched[[direction]]@features
	annotations = rna_patient_enrichment_go$enrich@annotation@annotation_features
	
	group_genes_metabolites = purrr::map(split_go, \(in_go){
		go_set = purrr::map(annotations[in_go], base::intersect, rna_sig)
		go_set = purrr::map(go_set, base::intersect, cor_genes)
		go_genes = unique(unlist(go_set))
		
		go_genes_cor = rna_within_cor |>
			dplyr::filter((s1 %in% go_genes) & (s2 %in% go_genes))
		gene_metabolite = sig_cor |>
			dplyr::filter(gene %in% go_genes)
		metabolite_cor = metabolites_within_cor |>
			dplyr::filter((s1 %in% unique(gene_metabolite$metabolite)) & (s2 %in% unique(gene_metabolite$metabolite)))
		
		list(go_set = go_set,
				 gene_cor = go_genes_cor,
				 metabolite_cor = metabolite_cor,
				 gene_metabolite = gene_metabolite)
	})
	
	group_genes_metabolites
}

annotate_lipids = function(lipidomics_keep)
{
	# tar_load(lipidomics_keep)
	lipid_data_df = rowData(lipidomics_keep) |> as.data.frame()
	fixed_lipids = fix_lipids(lipid_data_df$metabolite_name)
	annotated_lipids = rmf_annotate_lipids(fixed_lipids$Molecule)
	
	out_data = dplyr::left_join(fixed_lipids, annotated_lipids, by = "Molecule")
	
	lipids_2_annotation = dplyr::left_join(lipid_data_df[, c("feature_id", "metabolite_name")],
																				 out_data, by = c("metabolite_name" = "lipid"))
	
	lipids_2_annotation = lipids_2_annotation |>
		dplyr::filter(!is.na(not_matched) | !not_matched)
	
	lipids_2_annotation = lipids_2_annotation |>
		dplyr::mutate(class = class_stub,
									Class = NULL,
									total_length = total_cl,
									total_db = total_cs,
									chain1_length = l_1,
									chain2_length = l_2,
									chain3_length = l_3,
									chain1_db = s_1,
									chain2_db = s_2,
									chain3_db = s_3,
									class_total_length = paste0(class_stub, "_", total_length),
									class_total_db = paste0(class_stub, "_", total_db),
									class_chain1_length = paste0(class_stub, "_", chain1_length),
									class_chain2_length = paste0(class_stub, "_", chain2_length),
									class_chain3_length = paste0(class_stub, "_", chain3_length))
	
	annotation_cols = c("total_length", "total_db", 
											"chain1_length", "chain2_length", "chain3_length", "chain1_db", 
											"chain2_db", "chain3_db", "class_total_length", "class_total_db", 
											"class_chain1_length", "class_chain2_length", "class_chain3_length", 
											"class")
	
	annotation_list = purrr::map(annotation_cols, \(in_col){
		tmp_list = split(lipids_2_annotation$feature_id, lipids_2_annotation[[in_col]])
		names(tmp_list) = paste0(in_col, ":", names(tmp_list))
		tmp_list
	})
	
	feature_annotations = unlist(annotation_list, recursive = FALSE)
	na_annotations = grepl("_NA$", names(feature_annotations))
	feature_annotations = feature_annotations[!na_annotations]
	feature_annotations = purrr::map(feature_annotations, unique)
	
	lipid_annotation = categoryCompare2::annotation(feature_annotations, annotation_type = "lipids", feature_type = "lipids")
	
	lipid_annotation
}

choose_or = function(in_lipid)
{
	out_lipid = in_lipid |>
		stringr::str_split(pattern = "\\|") |>
		purrr::map_chr(\(in_split){
			if (length(in_split) > 1) {
				return(in_split[2])
			} else {
				return(in_split)
			}
		})
	out_lipid
}



fix_lipids = function(lipid_ids)
{
	# lipid_ids = "TG 48:3;O|TG 17:1_17:1_14:1;O"
	lipids_nona = lipid_ids[!is.na(lipid_ids)]
	
	lipid_ids2 = choose_or(lipids_nona) |>
		stringr::str_replace_all(" O-", "O ") |>
		stringr::str_replace_all("-O", "O") |>
		stringr::str_replace_all(" P-", "P ") |>
		stringr::str_replace_all("-N", "N") |>
		stringr::str_replace_all("Isomer [A|B|C|D]", "") |>
		stringr::str_replace_all("\\:\\:", ":")
	
	original_ids = data.frame(lipid = lipids_nona,
														pass1 = lipid_ids2)
	
	split_space = stringr::str_split(lipid_ids2, " ")
	names(split_space) = lipid_ids2
	has_o3 = grepl(";O3|;3O", names(split_space))
	
	split_o3 = split_space[has_o3]
	
	split_space = split_space[!has_o3]
	
	split_o3 = purrr::map(split_o3, \(in_o3){
		in_o3[1] = paste0(in_o3[1], "O3")
		in_o3[2] = gsub(";O3|;3O", "", in_o3[2])
		in_o3
	})
	
	has_o2 = grepl(";O2|;2O|\\(2OH\\)", names(split_space))
	split_o2 = split_space[has_o2]
	split_o2 = purrr::map(split_o2, \(in_o2){
		in_o2[1] = paste0(in_o2[1], "O2")
		in_o2[2] = gsub(";O2|;2O|;\\(2OH\\)|;O", "", in_o2[2])
		in_o2
	})
	
	split_space = split_space[!has_o2]
	
	has_o = grepl(";O", names(split_space))
	
	split_o = split_space[has_o]
	split_o = purrr::map(split_o, \(in_o){
		in_o[1] = paste0(in_o[1], "O")
		in_o[2] = gsub(";O", "", in_o[2])
		in_o
	})
	split_space = split_space[!has_o]
	split_all = c(split_space, split_o, split_o2, split_o3)
	
	split_all = purrr::map(split_all, \(in_split){
		in_split[1] = gsub("-", "", in_split[1])
		in_split
	})
	
	merged = purrr::map_chr(split_all, \(in_split){
		paste0(in_split, collapse = " ")
	})
	
	merged_df = data.frame(lipid_id = names(merged),
												 Molecule = merged)
	
	out_df = dplyr::left_join(original_ids, merged_df, by = c("pass1" = "lipid_id"),
														relationship = "many-to-many")
	out_df = out_df |>
		dplyr::select(-pass1)
	out_df
}

rmf_annotate_lipids = function(in_lipids)
{
	#in_lipids = fixed_lipids$Molecule
	annotated_lipids = lipidr::annotate_lipids(in_lipids)
	
	is_annotated = annotated_lipids |>
		dplyr::filter(!not_matched)
	
	has_multi = c("BMP", "Cer", "CerADS",
								"CerAS", "CerNDS", "CerO2", "CerO3",
								"DG", "FAHFA", "GalCer", "GlcCer",
								"Hex3CerO2", "HexCerNS", "HexCerO3",
								"MGDG", "PC", "PCO", "PE", "PECerO2",
								"PEO", "PEtOH", "PG", "PI",
								"PMeOH", "PS",
								"SHexCer", "SMO2", "TG", "TGO")
	is_annotated_should_multi = is_annotated |>
		dplyr::filter(class_stub %in% has_multi)
	
	is_annotated_whatever = is_annotated |>
		dplyr::filter(!(class_stub %in% has_multi))
	
	should_none_l1 = is_annotated_should_multi$l_1 == is_annotated_should_multi$total_cl
	is_annotated_should_multi[should_none_l1, "chain1"] = ""
	is_annotated_should_multi[should_none_l1, "l_1"] = NA
	is_annotated_should_multi[should_none_l1, "s_1"] = NA
	
	out_annotated = dplyr::bind_rows(is_annotated_should_multi,
																	 is_annotated_whatever)
	out_annotated
}

annotate_metabolite_type = function(metabolomics_feature_list)
{
	# tar_load(metabolomics_feature_list)
	split_type = split(metabolomics_feature_list$feature_id, metabolomics_feature_list$type)
	split_type = purrr::map(split_type, unique)
	annotation_type = categoryCompare2::annotation(split_type, annotation_type = "type",
																								 feature_type = "metabolite")
	annotation_type
}

run_binomial = function(metabolomics_de_patient_list,
												lipid_annotations)
{
	# tar_load(c(metabolomics_de_patient_list,
	# 				 lipid_annotations))
	pos_features = metabolomics_de_patient_list |>
		dplyr::filter(log2FoldChange > 0) |>
		dplyr::pull(feature_id)
	neg_features = metabolomics_de_patient_list |>
		dplyr::filter(log2FoldChange < 0) |>
		dplyr::pull(feature_id)
	
	binomial_res = binomial_feature_enrichment(
		new("binomial_features", positivefc = pos_features,
				negativefc = neg_features, annotation = lipid_annotations),
		p_adjust = "BH",
		min_features = 6
	)
	res_df = as.data.frame(binomial_res@statistics@statistic_data)
	res_df$id = binomial_res@statistics@annotation_id
	
	if (length(binomial_res@annotation@description) > 0) {
		res_df$description = binomial_res@annotation@description[res_df$id]
	}
	
	list(enrichment = binomial_res,
			 stats = res_df)
}

enrich_genes_correlated_lipids = function(rna_correlated_interesting_lipids,
																					ensembl_go)
{
	# tar_load(c(rna_correlated_interesting_lipids,
	# 					 ensembl_go))
	
	enrich_each = purrr::map(rna_correlated_interesting_lipids$groups, \(in_group){
		tmp_enrich = hypergeometric_feature_enrichment(
			new("hypergeom_features", significant = in_group$transcript,
					universe = rna_correlated_interesting_lipids$universe,
					annotation = ensembl_go),
			p_adjust = "BH", min_features = 2
		)
		tmp_df = as.data.frame(tmp_enrich@statistics@statistic_data)
		tmp_df$id = rownames(tmp_df)
		tmp_df$description = ensembl_go@description[tmp_df$id]
		list(enrich = tmp_enrich,
				 stats = tmp_df)
	})
	
}

}