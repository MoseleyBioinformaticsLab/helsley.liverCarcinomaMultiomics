create_dotplots = function(
  metabolite_collapsed_norm,
  rna_collapsed_norm,
  heatmap_correlations,
  color_scales,
  save_dir = "gene_compound_plots"
) {
  # tar_load(c(
  #   metabolite_collapsed_norm,
  #   rna_collapsed_norm,
  #   heatmap_correlations,
  #   color_scales
  # ))
  # save_dir = "gene_compound_plots"

  unlink(here::here(save_dir), recursive = TRUE)
  metabolite_data = rowData(metabolite_collapsed_norm)
  metabolite_counts = assays(metabolite_collapsed_norm)$counts

  rna_data = rowData(rna_collapsed_norm)
  rna_counts = assays(rna_collapsed_norm)$counts

  sample_info = colData(rna_collapsed_norm) |> tibble::as_tibble()
  use_samples = base::intersect(
    colnames(rna_counts),
    colnames(metabolite_counts)
  )

  all_plots = purrr::map_chr(
    seq_len(nrow(heatmap_correlations)),
    \(i_row) {
      r_loc = heatmap_correlations$s1[i_row]
      m_loc = heatmap_correlations$s2[i_row]

      plot_tibble = tibble::tibble(
        rna = rna_counts[r_loc, use_samples] + 1,
        metabolite = metabolite_counts[m_loc, use_samples] + 1,
        sample_id = use_samples
      )
      plot_tibble = dplyr::left_join(
        plot_tibble,
        sample_info[, c("sample_id", "treatment")],
        by = "sample_id"
      )
      r_id = rna_data[r_loc, ]$name
      m_id = metabolite_data[m_loc, ]$metabolite_id

      padj = heatmap_correlations$padjust[i_row]
      cor_val = heatmap_correlations$cor[i_row]

      filename = paste0(
        r_id,
        "_v_",
        janitor::make_clean_names(m_id),
        "_p_",
        formatC(padj, digits = 2, format = "f"),
        "_c_",
        formatC(cor_val, digits = 2, format = "f")
      )

      out_plot = plot_tibble |>
        ggplot(aes(x = log2(rna), y = log2(metabolite), color = treatment)) +
        geom_point(size = 2, show.legend = FALSE) +
        scale_color_manual(values = color_scales$treatment_normal_cancer) +
        labs(x = r_id, y = m_id)

      full_file = fs::path(here::here(save_dir), filename) |>
        paste0(".png")
      ragg::agg_png(full_file)
      print(out_plot)
      dev.off()
      filename
    }
  )

  heatmap_correlations = dplyr::left_join(
    heatmap_correlations,
    as.data.frame(rna_data[, c("feature_id", "name")]),
    by = c("s1" = "feature_id")
  )

  heatmap_correlations = dplyr::left_join(
    heatmap_correlations,
    as.data.frame(metabolite_data[, c("feature_id", "metabolite_id")]),
    by = c("s2" = "feature_id")
  )
  heatmap_correlations$file = all_plots

  return(heatmap_correlations)
}

filter_correlations = function(
  rna_compounds_matrix,
  rna_metabolites_all_spearman
) {
  rna_lipids = list(
    rna = rownames(rna_compounds_matrix$lipids),
    met = colnames(rna_compounds_matrix$lipids)
  )
  rna_met = list(
    rna = rownames(rna_compounds_matrix$compounds),
    met = colnames(rna_compounds_matrix$compounds)
  )

  rna_lipids_sub =
    lipids_correlations = rna_metabolites_all_spearman |>
      dplyr::filter((s1 %in% rna_lipids$rna) & (s2 %in% rna_lipids$met)) |>
      dplyr::mutate(heatmap_set = "lipids")

  met_correlations = rna_metabolites_all_spearman |>
    dplyr::filter((s1 %in% rna_met$rna) & (s2 %in% rna_met$met)) |>
    dplyr::mutate(heatmap_set = "metabolites")

  keep_correlations = dplyr::bind_rows(lipids_correlations, met_correlations)
  return(keep_correlations)
}

write_heatmap_file = function(heatmap_individual_plots, filename) {
  openxlsx::write.xlsx(
    heatmap_individual_plots,
    file = filename,
    overwrite = TRUE
  )
  return(filename)
}
