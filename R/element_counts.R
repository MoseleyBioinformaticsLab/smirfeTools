#' get IMF element counts
#'
#' Given an isotopic molecular formula (IMF), extract the counts for each element.
#'
#' @param imf_string the imf string to use
#' @return data.frame with element and count.
#'
#' @export
get_imf_elements <- function(imf_string){
  element_regex <- "[A-Z][a-z]?\\d*"
  numeric_regex <- "[0-9]+"

  extract_elements_imf <- stringr::str_locate_all(imf_string, element_regex)[[1]]
  extract_counts_imf <- stringr::str_locate_all(imf_string, numeric_regex)[[1]]

  imf_end_matches <- match(extract_counts_imf[,2],  extract_elements_imf[,2])
  extract_counts_imf <- extract_counts_imf[!is.na(imf_end_matches),]
  imf_end_matches <- imf_end_matches[!is.na(imf_end_matches)]
  extract_counts_imf <- extract_counts_imf[imf_end_matches, ]

  element_values_imf <- purrr::map_chr(seq(1, nrow(extract_elements_imf)),
                                       function(irow){
                                         substr(imf_string, extract_elements_imf[irow, 1], extract_counts_imf[irow, 1] - 1)
                                       })
  count_values_imf <- purrr::map_dbl(seq(1, nrow(extract_elements_imf)),
                                     function(irow){
                                       as.numeric(substr(imf_string, extract_counts_imf[irow, 1], extract_counts_imf[irow, 2]))
                                     })
  uniq_elements_imf <- unique(element_values_imf)
  uniq_counts_imf <- purrr::map_df(uniq_elements_imf, function(x){
    data.frame(element = x,
               count = sum(count_values_imf[element_values_imf %in% x]),
               stringsAsFactors = FALSE)
  })
  uniq_counts_imf
}


#' get EMF elements
#'
#' For an elemental molecular formula, extract the elements and their counts
#'
#' @param emf_string the EMF string to extract from.
#' @return data.frame with element and count
#' @export
get_emf_elements <- function(emf_string){
  element_regex <- "[A-Z][a-z]?\\d*"
  numeric_regex <- "[0-9]+"

  extract_elements_emf <- stringr::str_locate_all(emf_string, element_regex)[[1]]
  extract_counts_emf <- stringr::str_locate_all(emf_string, numeric_regex)[[1]]

  emf_end_matches <- match(extract_elements_emf[,2], extract_counts_emf[,2])

  element_values_emf <- purrr::map_chr(seq(1, nrow(extract_elements_emf)),
                                       function(irow){
                                         start <- extract_elements_emf[irow, 1]
                                         end <- emf_end_matches[irow]
                                         if (is.na(end)) {
                                           end <- extract_elements_emf[irow, 2]
                                         } else {
                                           end <- extract_counts_emf[emf_end_matches[irow], 1] - 1
                                         }
                                         substr(emf_string, start, end)
                                       })

  element_counts_emf <- purrr::map_dbl(seq(1, nrow(extract_elements_emf)),
                                       function(irow){
                                         if (is.na(emf_end_matches[irow])) {
                                           count <- 1
                                         } else {
                                           end_match <- emf_end_matches[irow]
                                           count <- as.numeric(substr(emf_string, extract_counts_emf[end_match, 1], extract_counts_emf[end_match, 2]))
                                         }
                                         count
                                       })



  uniq_elements_emf <- unique(element_values_emf)
  uniq_counts_emf <- purrr::map_df(uniq_elements_emf, function(x){
    data.frame(element = x,
               count = sum(element_counts_emf[element_values_emf %in% x]),
               stringsAsFactors = FALSE)
  })
  uniq_counts_emf <- uniq_counts_emf[order(uniq_counts_emf$element), ]
  rownames(uniq_counts_emf) <- NULL
  uniq_counts_emf
}

#' Compare IMF and EMF Counts
#'
#' Given an IMF counts `data.frame` and an EMF counts `data.frame`, compare them
#' to make sure that they contain identical counts of each element.
#'
#' @param imf_counts data.frame from [get_imf_elements]
#' @param emf_counts data.frame from [get_emf_elements]
#'
#' @export
#' @return logical
check_imf_emf_elemental_counts <- function(imf_counts, emf_counts){
  joined_counts <- dplyr::full_join(imf_counts, emf_counts, by = "element", suffix = c(".imf", ".emf"))
  if ((sum(is.na(joined_counts$element)) == 0) && identical(joined_counts$count.imf, joined_counts$count.emf)) {
    TRUE
  } else {
    FALSE
  }
}