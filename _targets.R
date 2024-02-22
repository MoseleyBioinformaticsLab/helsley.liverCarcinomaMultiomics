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
		(\(x){setup_metabolomics(x, sample_info)})(),
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
		(\(x){setup_metabolomics(x, sample_info)})(),
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
		(\(x){setup_metabolomics(x, sample_info)})(),
	primary_metabolism_norm = median_normalization(primary_metabolism),
	primary_metabolism_keep = keep_presence(primary_metabolism_norm),
	pr_n_qcqa = c(get_n_features(primary_metabolism_norm), get_n_features(primary_metabolism_keep)),
	
	primary_metabolism_cor_pca = sample_correlations_pca(primary_metabolism_keep),
	primary_metabolism_qcqa = create_qcqa_plots(primary_metabolism_cor_pca,
																							color_scales),
	
	## Differential Analysis --------
	### Determine Outliers --------
	rna_outliers = determine_outliers(rna_keep),
	bioamines_outliers = determine_outliers(bioamines_keep),
	lipidomics_outliers = determine_outliers(lipidomics_keep),
	pr_outliers = determine_outliers(primary_metabolism_keep),
	
	### Collapse Replicates -----
	### Note that this function also *removes* the previously identified outliers before
	### collapsing the replicates.
	rna_collapsed = collapse_deseq_replicates(rna_outliers),
	bioamines_collapsed = collapse_metabolomics_replicates(bioamines_outliers),
	lipidomics_collapsed = collapse_metabolomics_replicates(lipidomics_outliers),
	pr_collapsed = collapse_metabolomics_replicates(pr_outliers),
	
	### Unpaired --------
	rna_de_treatment = calculate_deseq_stats(rna_collapsed,
																					 which = "treatment"),
	bioamines_de_treatment = calculate_metabolomics_stats(bioamines_collapsed,
																												which = "unpaired"),
	lipidomics_de_treatment = calculate_metabolomics_stats(lipidomics_collapsed,
																												 which = "unpaired"),
	pr_de_treatment = calculate_metabolomics_stats(pr_collapsed,
																								 which = "unpaired"),
	### Paired --------
	rna_paired = filter_to_pairs(rna_collapsed),
	rna_de_patient = calculate_deseq_stats(rna_paired,
																				 which = "patient"),
	
	bioamines_paired = filter_to_pairs(bioamines_collapsed),
	
	## documents -----------
	tar_quarto(qcqa, "docs/qcqa.qmd")
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
