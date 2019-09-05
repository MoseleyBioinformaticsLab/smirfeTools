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