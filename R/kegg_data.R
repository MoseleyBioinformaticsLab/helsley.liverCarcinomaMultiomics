get_kegg_compound_data = function(){
  compound_2_pathway = link_df("compound", "pathway") |>
    dplyr::mutate(pathway = gsub("path:", "", pathway),
    							compound = gsub("cpd:", "", compound))
  
  gene_2_pathway = link_df("hsa", "pathway") |>
  	dplyr::mutate(pathway = gsub("hsa", "map", gsub("path:", "", pathway)),
  								gene = gsub("hsa:", "", hsa))

  compound_2_enzyme = link_df("compound", "enzyme")

  enzyme_2_gene = link_df("enzyme", "hsa")

  pathway_id = list_df("pathway")

  pathway_human_id = list_df("pathway", "hsa")

  compound_2_pathway = compound_2_pathway |>
    dplyr::filter(!(grepl("map01100", pathway)))

  compound_2_enzyme = compound_2_enzyme |>
    dplyr::filter(enzyme %in% enzyme_2_gene$enzyme)

  compound_id = list_df("compound")

  enzyme_id = list_df("enzyme")

  # compound_id = compound_id %>%
  #   dplyr::filter(compound %in% compound_2_enzyme$compound)
  # enzyme_id = enzyme_id %>%
  #   dplyr::filter(enzyme %in% compound_2_enzyme$enzyme)
  list(compound_2_pathway = compound_2_pathway,
       compound_2_enzyme = compound_2_enzyme,
  		 gene_2_pathway = gene_2_pathway,
       enzyme_2_gene = enzyme_2_gene,
       pathway_human_id = pathway_human_id,
       pathway_id = pathway_id,
       compound_id = compound_id,
       enzyme_id = enzyme_id)
}

link_df = function(subject, target){
  subject_2_target_list = KEGGREST::keggLink(subject, target)

  sub_2_target_df = data.frame(subject = subject_2_target_list,
                               target = names(subject_2_target_list))
  names(sub_2_target_df) = c(subject, target)
  rownames(sub_2_target_df) = NULL
  sub_2_target_df
}

list_df = function(subject, target = NULL){
  if (!is.null(target)) {
    subject_2_target_list = KEGGREST::keggList(subject, target)
    sub_2_target_df = data.frame(subject = subject_2_target_list,
                                 target = names(subject_2_target_list))
    names(sub_2_target_df) = c(subject, target)
  } else {
    subject_2_target_list = KEGGREST::keggList(subject)
    sub_2_target_df = data.frame(subject = subject_2_target_list,
                                 target = names(subject_2_target_list))
    names(sub_2_target_df) = c("description", subject)
  }
  rownames(sub_2_target_df) = NULL
  sub_2_target_df
}

create_kegg_annotations = function(feature_metadata,
                                   kegg_data)
{
  # tar_load(feature_metadata)
  # tar_load(kegg_data)
  human_maps = kegg_data$pathway_human_id |>
    dplyr::mutate(pathway_description = gsub("-.*", "", pathway),
                  pathway = gsub("hsa", "map", hsa))
  compounds = kegg_data$compound_2_pathway |>
    dplyr::mutate(compound = gsub("cpd:", "", compound)) |>
    dplyr::filter(pathway %in% human_maps$pathway)
  split_maps = split(compounds$compound, compounds$pathway)
  split_maps = purrr::map(split_maps, unique)
  pathway_descriptions = human_maps$pathway_description
  names(pathway_descriptions) = human_maps$pathway
  pathway_descriptions = pathway_descriptions[names(split_maps)]
  kegg_annotation = categoryCompare2::annotation(split_maps,
                                                 annotation_type = "kegg",
                                                 description = pathway_descriptions,
                                                 feature_type = "metabolite")
  kegg_annotation
}

create_kegg_annotations_genes = function(kegg_data,
																				 ensembl_entrez)
{
	# tar_load(kegg_data)
	# tar_load(ensembl_entrez)
	human_maps = kegg_data$pathway_human_id |>
		dplyr::mutate(pathway_description = gsub("-.*", "", pathway),
									pathway = gsub("hsa", "map", hsa))
	genes = kegg_data$gene_2_pathway |>
		dplyr::filter(pathway %in% human_maps$pathway)
	genes = dplyr::left_join(genes, ensembl_entrez |>
													 	dplyr::mutate(gene = as.character(entrezgene_id)), by = "gene",
													 relationship = "many-to-many")
	split_maps = split(genes$ensembl_gene_id, genes$pathway)
	split_maps = purrr::map(split_maps, unique)
	pathway_descriptions = human_maps$pathway_description
	names(pathway_descriptions) = human_maps$pathway
	pathway_descriptions = pathway_descriptions[names(split_maps)]
	kegg_annotation = categoryCompare2::annotation(split_maps,
																								 annotation_type = "kegg",
																								 description = pathway_descriptions,
																								 feature_type = "ENSEMBL")
	kegg_annotation
}

create_kegg_annotations_compounds = function(kegg_data,
																						 feature_to_kegg_chebi)
{
	# tar_load(kegg_data)
	# tar_load(feature_to_kegg_chebi)
	human_maps = kegg_data$pathway_human_id |>
		dplyr::mutate(pathway_description = gsub("-.*", "", pathway),
									pathway = gsub("hsa", "map", hsa))
	compounds = kegg_data$compound_2_pathway |>
		dplyr::filter(pathway %in% human_maps$pathway)
	compounds = dplyr::left_join(compounds, feature_to_kegg_chebi, by = c("compound" = "kegg"), 
															 relationship = "many-to-many")
	split_maps = split(compounds$feature_id, compounds$pathway)
	split_maps = purrr::map(split_maps, unique)
	pathway_descriptions = human_maps$pathway_description
	names(pathway_descriptions) = human_maps$pathway
	pathway_descriptions = pathway_descriptions[names(split_maps)]
	kegg_annotation = categoryCompare2::annotation(split_maps,
																								 annotation_type = "kegg",
																								 description = pathway_descriptions,
																								 feature_type = "feature_id")
	kegg_annotation
}
