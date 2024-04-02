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