create_dotplots = function(
  rna_compounds_matrix,
  metabolite_collapsed_norm,
  rna_collapsed_norm,
  rna_metabolites_all_spearman,
  color_scales
) {
  # tar_load(c(
  #   rna_compounds_matrix,
  #   metabolite_collapsed_norm,
  #   rna_collapsed_norm,
  #   rna_metabolites_all_spearman,
  #   color_scales
  # ))

  return(NULL)
}

filter_correlations = function(
  rna_compounds_matrix,
  rna_metabolites_all_spearman
) {
  all_ens = c(
    rownames(rna_compounds_matrix$compounds),
    rownames(rna_compounds_matrix$lipids)
  )
  all_comp = c(
    colnames(rna_compounds_matrix$compounds),
    colnames(rna_compounds_matrix$lipids)
  )

  keep_correlations = rna_metabolites_all_spearman |>
    dplyr::filter((s1 %in% all_ens) & (s2 %in% all_comp))

  return(keep_correlations)
}
