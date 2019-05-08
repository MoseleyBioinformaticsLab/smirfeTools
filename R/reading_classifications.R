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

  zero_frame = data.frame(Categories = NA, Classes = NA, isotopologue_EMF = "NA",
                          stringsAsFactors = FALSE)
  zero_frame = zero_frame[0, ]

  has_category = purrr::map_lgl(emf_classes, ~ length(.x$Categories) > 0)
  process_classes = emf_classes[has_category]
  class_data <- purrr::map2_dfr(process_classes, names(process_classes), function(in_list, in_emf){
    #message(in_emf)
    n_category = length(in_list$Categories)
    n_classes = length(in_list$Classes)

    if (n_classes == 0) {
      in_list$Classes = NA
    }

    if (((n_category > 0) && (n_category == n_classes)) || (n_classes == 0)) {
      tmp_frame = as.data.frame(in_list, stringsAsFactors = FALSE)
      tmp_frame$isotopologue_EMF = in_emf
    } else if (n_category > 0) {
      tmp_frame = as.data.frame(in_list, stringsAsFactors = FALSE)
      tmp_frame$isotopologue_EMF = in_emf
      #message(in_emf)
    } else {
      tmp_frame = zero_frame
    }

    tmp_frame
  })

  isnt_not_or_na = !((class_data$Categories %in% "not_lipid") | (is.na(class_data$Classes)))

  tmp_class = class_data[!isnt_not_or_na, ]

  has_class = class_data[isnt_not_or_na, ]
  extract_pattern = "\\[.*\\]"
  check_class = stringr::str_extract(has_class$Categories, extract_pattern) %>% gsub("\\[|\\]", "", .)

  has_class2 = purrr::map_df(seq(1, nrow(has_class)), function(in_row){
    tmp_row = has_class[in_row, ]
    if (!(grepl(check_class[in_row], tmp_row[1, "Classes"]))) {
      tmp_row[1, "Classes"] = "NA"
    }
    tmp_row
  })
  out_classes = rbind(tmp_class, has_class2)
  out_classes
}
