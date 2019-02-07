#' Import EMF Classifications
#'
#' Allows the reading in of EMF classification results from a JSON file.
#'
#' @param classification_file the file to read from
#' @param fix_json properly format the JSON so the JSON parser can read it?
#' @param rm_na_formula remove those entries that have a formula entry of `NA`??
#'
#' @return data.frame
#' @export
import_emf_classifications = function(classification_file, fix_json = TRUE, rm_na_formula = TRUE){
  if (fix_json) {
    emf_json = gsub("None", "null",
                    gsub("'", '"',
                         scan(classification_file, what = character(), sep = "\n", quiet = TRUE)
                    )
    )
    emf_classes = purrr::map(emf_json, jsonlite::fromJSON)
  } else {
    emf_classes = jsonlite::fromJSON(classification_file)
  }

  class_data_names = unique(unlist(purrr::map(emf_classes, names)))

  class_list = as.list(rep(NA, length(class_data_names)))
  names(class_list) = class_data_names
  class_tmp_df = as.data.frame(class_list)

  class_data <- purrr::map_df(emf_classes, function(in_list){
    #message(in_list$SortedFormula)
    length_data <- purrr::map_int(in_list, length)
    max_length = max(length_data)
    in_list = in_list[length_data > 0]

    in_list = purrr::map(in_list, function(x){
      if (length(x) == 1) {
        rep(x, max_length)
      } else {
        sort(c(x, rep(NA, max_length - length(x))), na.last = TRUE)
      }
    })
    tmp_df = purrr::map_df(seq(1, max_length), ~ class_tmp_df)
    tmp_df[, names(in_list)] <- in_list
    #print(warnings())
    tmp_df
  })
  if (rm_na_formula) {
    class_data = dplyr::filter(class_data, !(is.na(Formula)))
  }

  class_data
}
