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
	
	rna_norm = DESeq2::counts(rna_dds, normalized = TRUE) |>
		matrix_2_long() |>
		floor_values(),
	rna_keep = keep_presence(rna_norm,
													 sample_info),
	rna_n_qcqa = c(nrow(rna_dds), nrow(rna_keep)),
	
	rna_ratios = calculate_patient_ratios(rna_keep,
																				sample_info_no11),
	
	rna_cor_pca = sample_correlations_pca(rna_keep, sample_info),
	rna_qcqa = create_qcqa_plots(rna_cor_pca,
															 sample_info,
															 color_scales),
	
	rna_ratios_cor_pca = sample_correlations_pca(rna_ratios, patient_info),
	rna_ratios_qcqa = create_qcqa_plots(rna_ratios_cor_pca,
																			patient_info,
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
		split_intensities_feature_metadata(),
	bioamines_norm = median_normalization(bioamines$abundance),
	bioamines_keep = keep_presence(bioamines_norm,
																 sample_info),
	
	bioamines_n_qcqa = c(nrow(bioamines_norm), nrow(bioamines_keep)),
	bioamines_ratios = calculate_patient_ratios(bioamines_keep, sample_info_no11),
	
	bioamines_cor_pca = sample_correlations_pca(bioamines_keep, sample_info),
	bioamines_qcqa = create_qcqa_plots(bioamines_cor_pca,
																		 sample_info,
																		 color_scales),
	
	bioamines_ratios_cor_pca = sample_correlations_pca(bioamines_ratios,
																										 patient_info),
	bioamines_ratios_qcqa = create_qcqa_plots(bioamines_ratios_cor_pca,
																						patient_info,
																						color_scales),
	
	### Lipidomics -----
	tar_target(lipidomics_file,
						 "raw_data/metabolomics_lipidomics.xlsx",
						 format = "file"),
	lipidomics = suppressMessages(readxl::read_excel(lipidomics_file,
																									 sheet = "Data",
																									 skip = 8)) |>
		janitor::clean_names()  |>
		rename_experimental_samples() |>
		split_intensities_feature_metadata(),
	lipidomics_norm = median_normalization(lipidomics$abundance),
	lipidomics_keep = keep_presence(lipidomics_norm,
																	sample_info),
	lipidomics_n_qcqa = c(nrow(lipidomics_norm), nrow(lipidomics_keep)),
	lipidomics_ratios = calculate_patient_ratios(lipidomics_keep,
																							 sample_info_no11),
	
	lipidomics_cor_pca = sample_correlations_pca(lipidomics_keep, sample_info),
	lipidomics_qcqa = create_qcqa_plots(lipidomics_cor_pca,
																			sample_info,
																			color_scales),
	lipidomics_ratios_cor_pca = sample_correlations_pca(lipidomics_ratios,
																											patient_info),
	lipidomics_ratios_qcqa = create_qcqa_plots(lipidomics_ratios_cor_pca,
																						 patient_info,
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
		split_intensities_feature_metadata(),
	primary_metabolism_norm = median_normalization(primary_metabolism$abundance),
	primary_metabolism_keep = keep_presence(primary_metabolism_norm,
																					sample_info),
	pr_n_qcqa = c(nrow(primary_metabolism_norm), nrow(primary_metabolism_keep)),
	primary_metabolism_ratios = calculate_patient_ratios(primary_metabolism_keep,
																											 sample_info_no11),
	
	primary_metabolism_cor_pca = sample_correlations_pca(primary_metabolism_keep,
																											 sample_info),
	primary_metabolism_qcqa = create_qcqa_plots(primary_metabolism_cor_pca,
																							sample_info,
																							color_scales),
	primary_metabolism_ratios_cor_pca = sample_correlations_pca(primary_metabolism_ratios,
																														 patient_info),
	primary_metabolism_ratios_qcqa = create_qcqa_plots(primary_metabolism_ratios_cor_pca,
																										 patient_info,
																										 color_scales),
	
	## Differential Analysis --------
	
	### Unpaired --------
	#### RNA -------
	# rna_dds_treatment = prep_rna_data_treatment(rna_data,
	# 																						sample_info_no11),
	# rna_dds_keep_treatment = keep_presence_dds(rna_dds_treatment, 0.75),
	# rna_deseq2_treatment = DESeq2::DESeq(rna_dds_keep_treatment),
	# rna_results_treatment = DESeq2::results(rna_deseq2_treatment, contrast = c("treatment", "cancerous", "normal_adjacent"), tidy = TRUE),
	## documents -----------
	tar_quarto(qcqa, "docs/qcqa.qmd")
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
