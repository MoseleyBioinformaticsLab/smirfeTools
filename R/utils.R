na_dataframe = function(template_dataframe){
  df_types = purrr::map_chr(template_dataframe, class)
  names(df_types) = names(template_dataframe)

  na_df = template_dataframe[1, , drop = FALSE]

  for (itype in names(df_types)) {
    na_df[[itype]] = as(NA, df_types[itype])
  }
  rownames(na_df) = NULL
  na_df
}

check_ngemf = function(in_emf, min_emfs = 2){
  if (!is.null(in_emf)) {
    n_gemf = length(unique(in_emf$grouped_EMF))
    if (n_gemf >= min_emfs) {
      return(in_emf)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}
