## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
tar_source("R")

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

	tar_target(sample_list_file,
						 "raw_data/sample_list.xlsx",
						 format = "file"),
	
	sample_info = suppressMessages(readxl::read_excel(sample_list_file)) |>
		janitor::clean_names() |>
		dplyr::mutate(sample_id = paste0("s", sample_label),
									treatment = tolower(gsub(" ", "_", treatment)),
									patient = tolower(gsub(" ", "_", patient))),
	
	## RNA -------
	tar_target(rna_file,
						 "raw_data/transcriptomics_counts.txt",
						 format = "file"),
	rna_data = suppressMessages(readr::read_tsv(rna_file)) |>
		janitor::clean_names(),
	
	rna_dds = prep_rna_data(rna_data[, seq_len(16)], sample_info),
	rna_norm = DESeq2::counts(rna_dds, normalized = TRUE) |>
		matrix_2_long(),
	rna_keep = keep_presence(rna_norm,
													 sample_info),
	
	rna_cor_pca = sample_correlations(rna_keep),
	
	
	## Biogenic Amines -----
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
	bioamines_cor_pca = sample_correlations(bioamines_keep),
	
	## Lipidomics -----
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
	lipidomics_cor_pca = sample_correlations(lipidomics_keep),
	
	## Primary Metabolism ----
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
	primary_metabolism_cor_pca = sample_correlations(primary_metabolism_keep)
	
	
# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style

)
