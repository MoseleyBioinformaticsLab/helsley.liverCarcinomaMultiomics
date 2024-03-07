## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
tar_source("R")

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

	color_scales = sample_colors(),
	
	tar_target(sample_list_file,
						 "raw_data/sample_list.xlsx",
						 format = "file"),
	
	sample_info = suppressMessages(readxl::read_excel(sample_list_file)) |>
		janitor::clean_names() |>
		dplyr::mutate(sample_id = paste0("s", sample_label),
									treatment = tolower(gsub(" ", "_", treatment)),
									patient = tolower(gsub(" ", "_", patient))),
	
	sample_info_no11 = sample_info |>
		dplyr::filter(!(sample_id %in% "s11")),
	patient_info = sample_info_no11 |>
		dplyr::transmute(sample_id = patient,
									treatment = "none") |>
		dplyr::distinct(),
	## QC-QA ------
	### RNA -------
	tar_target(rna_file,
						 "raw_data/transcriptomics_counts.txt",
						 format = "file"),
	rna_data = suppressMessages(readr::read_tsv(rna_file)) |>
		janitor::clean_names(),
	
	rna_dds = prep_rna_data_treatment(rna_data = rna_data, sample_metadata = sample_info,
																	count_locs = seq(2, 16),
																	feature_metadata = seq(17, 25)),
	
	novagene_rna_stats = suppressMessages(readr::read_tsv("raw_data/rnaseq_control_vs_tumor.txt")) |>
		janitor::clean_names() |>
		dplyr::select(gene_id, log2fold_change, pvalue, padj),
	
	rna_keep = keep_presence_dds(rna_dds),
	rna_n_qcqa = c(get_n_features(rna_dds), get_n_features(rna_keep)),
	
	rna_cor_pca = sample_correlations_pca(rna_keep),
	rna_qcqa = create_qcqa_plots(rna_cor_pca,
															 color_scales),
	
	### Biogenic Amines -----
	tar_target(bioamines_file,
						 "raw_data/metabolomics_biogenic_amines.xlsx",
						 format = "file"),
	bioamines = suppressMessages(readxl::read_excel(bioamines_file,
																									sheet = "Data",
																									skip = 9)) |>
		janitor::clean_names() |>
		rename_experimental_samples() |>
		(\(x){setup_metabolomics(x, sample_info, "bioamine")})(),
	bioamines_norm = median_normalization(bioamines),
	bioamines_keep = keep_presence(bioamines_norm),
	
	bioamines_n_qcqa = c(get_n_features(bioamines_norm), get_n_features(bioamines_keep)),
	
	bioamines_cor_pca = sample_correlations_pca(bioamines_keep),
	bioamines_qcqa = create_qcqa_plots(bioamines_cor_pca,
																		 color_scales),
	
	### Lipidomics -----
	tar_target(lipidomics_file,
						 "raw_data/metabolomics_lipidomics.xlsx",
						 format = "file"),
	lipidomics = suppressMessages(readxl::read_excel(lipidomics_file,
																									 sheet = "Data",
																									 skip = 8)) |>
		janitor::clean_names() |>
		rename_experimental_samples() |>
		(\(x){setup_metabolomics(x, sample_info, "lipids")})(),
	lipidomics_norm = median_normalization(lipidomics),
	lipidomics_keep = keep_presence(lipidomics_norm),
	lipidomics_n_qcqa = c(get_n_features(lipidomics_norm), get_n_features(lipidomics_keep)),
	
	lipidomics_cor_pca = sample_correlations_pca(lipidomics_keep),
	lipidomics_qcqa = create_qcqa_plots(lipidomics_cor_pca,
																			color_scales),
	
	### Primary Metabolism ----
	tar_target(primary_metabolism_file,
						 "raw_data/metabolomics_primary_metabolism.xlsx",
						 format = "file"),
	primary_metabolism = suppressMessages(readxl::read_excel(primary_metabolism_file,
																													 sheet = "data",
																													 skip = 8)) |>
		janitor::clean_names() |>
		rename_experimental_samples() |>
		(\(x){setup_metabolomics(x, sample_info, "pm")})(),
	primary_metabolism_norm = median_normalization(primary_metabolism),
	primary_metabolism_keep = keep_presence(primary_metabolism_norm),
	pm_n_qcqa = c(get_n_features(primary_metabolism_norm), get_n_features(primary_metabolism_keep)),
	
	primary_metabolism_cor_pca = sample_correlations_pca(primary_metabolism_keep),
	primary_metabolism_qcqa = create_qcqa_plots(primary_metabolism_cor_pca,
																							color_scales),
	
	## Differential Analysis --------
	### Determine Outliers --------
	rna_outliers = determine_outliers(rna_keep),
	bioamines_outliers = determine_outliers(bioamines_keep),
	lipidomics_outliers = determine_outliers(lipidomics_keep),
	pm_outliers = determine_outliers(primary_metabolism_keep),
	
	### Collapse Replicates -----
	### Note that this function also *removes* the previously identified outliers before
	### collapsing the replicates.
	rna_collapsed = collapse_deseq_replicates(rna_outliers),
	bioamines_collapsed = collapse_metabolomics_replicates(bioamines_outliers),
	lipidomics_collapsed = collapse_metabolomics_replicates(lipidomics_outliers),
	pm_collapsed = collapse_metabolomics_replicates(pm_outliers),
	
	### Unpaired --------
	rna_de_treatment = calculate_deseq_stats(rna_collapsed,
																					 which = "treatment"),
	bioamines_de_treatment = calculate_metabolomics_stats(bioamines_collapsed),
	lipidomics_de_treatment = calculate_metabolomics_stats(lipidomics_collapsed),
	pm_de_treatment = calculate_metabolomics_stats(pm_collapsed),
	### Paired --------
	rna_paired = filter_to_pairs(rna_collapsed),
	rna_de_patient = calculate_deseq_stats(rna_paired,
																				 which = "patient"),
	
	bioamines_paired = filter_to_pairs(bioamines_collapsed),
	bioamines_de_patient = calculate_metabolomics_stats(bioamines_paired, paired = "patient"),
	lipidomics_paired = filter_to_pairs(lipidomics_collapsed),
	lipidomics_de_patient = calculate_metabolomics_stats(lipidomics_paired, paired = "patient"),
	pm_paired = filter_to_pairs(pm_collapsed),
	pm_de_patient = calculate_metabolomics_stats(pm_paired,
																							 paired = "patient"),
	
	### Paired Samples, but unpaired stats --------
	rna_de_patient_unpaired = calculate_deseq_stats(rna_paired,
																				 which = "treatment"),
	
	bioamines_de_patient_unpaired = calculate_metabolomics_stats(bioamines_paired, paired = "treatment"),
	lipidomics_de_patient_unpaired = calculate_metabolomics_stats(lipidomics_paired, paired = "treatment"),
	pm_de_patient_unpaired = calculate_metabolomics_stats(pm_paired,
																							 paired = "treatment"),

	## Annotations --------
	### Genes -------
	ensembl_uniprot = get_ensembl_uniprot(version = "111"),
	tar_target(go_file,
						 "data/ancestors.json",
						 format = "file"),
	tar_target(namespace_file,
						 "data/namespace.json",
						 format = "file"),
	ensembl_go = create_go_annotations(ensembl_uniprot,
																		 go_file,
																		 namespace_file),
	tar_target(reactome_gene_file,
						 "data/Ensembl2Reactome_All_Levels.txt",
						 format = "file"),
	ensembl_reactome = create_reactome_gene_annotations(reactome_gene_file),
	
	### Metabolites -----
	tar_target(chebi_reactome_file,
						 "data/ChEBI2Reactome_All_Levels.txt",
						 format = "file"),
	chebi_reactome = create_reactome_chebi_annotations(chebi_reactome_file),
	tar_target(chebi_inchi_file,
						 "data/chebiId_inchi.tsv.gz",
						 format = "file"),
	inchi_df = write_inchi(chebi_inchi_file),
	inchikey_hash = get_inchikey_hash(inchi_df,
																		"data/chebi"),
	chebi_inchikey = map_chebi_to_inchikey(inchi_df,
																				 inchi_dir = "data/chebi",
																				 inchikey_hash),
	
	rna_treatment_enrichment_go = run_enrichment(rna_de_treatment, ensembl_go),
	rna_patient_enrichment_go = run_enrichment(rna_de_patient, ensembl_go),
	
	rna_treatment_enrichment_reactome = run_enrichment(rna_de_treatment, ensembl_reactome),
	rna_patient_enrichment_reactome = run_enrichment(rna_de_patient, ensembl_reactome),
	
	metabolomics_de_treatment_list = list(bioamines = bioamines_de_treatment,
																				lipidomics = lipidomics_de_treatment,
																				pm = pm_de_treatment),
	metabolomics_de_treatment = merge_and_map_metabolomics(metabolomics_de_treatment_list,
																												 chebi_inchikey),
	
	metabolomics_treatment_enrichment_reactome = run_enrichment(metabolomics_de_treatment,
																															chebi_reactome),
	
	
	metabolomics_de_patient_list = list(bioamines = bioamines_de_patient,
																				lipidomics = lipidomics_de_patient,
																				pm = pm_de_patient),
	metabolomics_de_patient = merge_and_map_metabolomics(metabolomics_de_patient_list,
																												 chebi_inchikey),
	
	metabolomics_patient_enrichment_reactome = run_enrichment(metabolomics_de_patient,
																															chebi_reactome),
	
	
	## Comparing Statistics Between Treatment and Paired ---------
	metabolomics_de_compare = compare_treatment_patient(metabolomics_de_treatment,
																											metabolomics_de_patient),
	
	rna_de_compare = compare_treatment_patient(rna_de_treatment,
																											rna_de_patient),
	
	metabolomics_de_patient_unpaired_list = list(bioamines = bioamines_de_patient_unpaired,
																			lipidomics = lipidomics_de_patient_unpaired,
																			pm = pm_de_patient_unpaired),
	metabolomics_de_patient_unpaired = merge_and_map_metabolomics(metabolomics_de_patient_unpaired_list,
																											 chebi_inchikey),
	rna_de_unpaired_compare = compare_treatment_patient(
		rna_de_patient,
		rna_de_patient_unpaired
	),
	
	metabolomics_de_unpaired_compare = compare_treatment_patient(
		metabolomics_de_patient,
		metabolomics_de_patient_unpaired
	),
	
	## documents -----------
	tar_quarto(qcqa, "docs/qcqa.qmd"),
	
	tar_quarto(de_comparisons, "docs/de_comparisons.qmd"),
	tar_quarto(differential_analysis, "docs/differential_analysis.qmd")
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
