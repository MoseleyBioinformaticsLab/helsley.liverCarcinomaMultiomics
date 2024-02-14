rename_experimental_samples = function(metabolomics_dataset)
{
	# metabolomics_dataset = tar_read(bioamines)
	names(metabolomics_dataset) = gsub("^x", "s", names(metabolomics_dataset))
	metabolomics_dataset
}