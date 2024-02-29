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


run_enrichment = function(de_values,
													annotation_obj,
													padj_cutoff = 0.01)
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
	enrich_comb = combine_enrichments(all = enrich_all,
																		pos = enrich_pos,
																		neg = enrich_neg)
	stats_list = list(all = stats_all,
										pos = stats_pos,
										neg = stats_neg)
	return(list(enrich = enrich_comb,
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
