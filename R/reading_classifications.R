#' Import EMF Classifications
#'
#' Allows the reading in of EMF classification results from a JSON file.
#'
#' @param classification_file the file to read from
#' @param fix_json properly format the JSON so the JSON parser can read it?
#'
#' @return data.frame
#' @export
import_emf_classifications = function(classification_file, fix_json = TRUE){
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
    length_data <- purrr::map_int(in_list, length)
    tmp_df = class_tmp_df
    if ((length_data["Categories"] != 0) || (length_data["PredictedCategories"] != 0)) {
      in_list <- in_list[length_data > 0]
      tmp_df[1, names(in_list)] = in_list
    }

    tmp_df
  })
  class_data = dplyr::filter(class_data, !(is.na(Formula)))
  class_data
}
