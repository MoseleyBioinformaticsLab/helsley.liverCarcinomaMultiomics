parse_kegg_compound = function(compound_file)
{
  # compound_file = "data/kegg_compounds/C00430.txt"
  compound_data = scan(file = compound_file, what = character(), sep = "\n",
                       quiet = TRUE)
  formula_or_comment = suppressWarnings(min(which(grepl("FORMULA|COMMENT|REACTION|REMARK|SEQUENCE", compound_data))))
  if (is.infinite(formula_or_comment)) {
    return(NULL)
  }

  info_loc = seq(1, formula_or_comment - 1)
  info_data = compound_data[info_loc]

  split_info = strsplit(info_data, "  |;")
  more_info = purrr::map(split_info, \(x){
    out_x = x[nchar(x) > 0]
    trimws(out_x)
  })
  n_info = length(more_info)
  first_synonym = more_info[[2]][2]
  if (length(more_info) > 2) {
    syn_locs = seq(3, n_info)
    more_synonym = purrr::map_chr(more_info[syn_locs], \(x){x})
    all_synonym = c(first_synonym, more_synonym)
  } else {
    all_synonym = first_synonym
  }
  info_out = tibble::tibble(kegg = more_info[[1]][2],
                            synonym = tolower(all_synonym))
  return(info_out)
}

get_all_kegg_compound = function()
{
  kegg_files = dir("data/kegg_compounds", pattern = "txt$", full.names = TRUE)
  kegg_compound_synonym = purrr::map(kegg_files, parse_kegg_compound) |>
    dplyr::bind_rows()
  kegg_compound_synonym
}
