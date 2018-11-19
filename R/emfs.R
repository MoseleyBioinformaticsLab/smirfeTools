#' EMFs to Peaks
#'
#' If we want to work with EMFs, we should first determine what are `unique` EMFs
#' based on the mapping between peaks and EMFs. This function adds a new column
#' to the data.frame, `grouped_EMF`
#'
#' @param peak_info data.frame mapping peaks to EMFs
#' @param peak_var which column has the peak information
#' @param emf_var which column has the EMF information
#' @param grouping_var what should the grouping variable be called
#'
#'
#' @return data.frame
#' @export
#'
#' @examples
group_emfs_by_peaks = function(peak_info, peak_var = "peak", emf_var = "complete_EMF",
                               grouping_var = "grouped_EMF"){
  split_emfs = split(peak_info[, peak_var], peak_info[, emf_var])

  emf_list = vector("list", length(split_emfs))

  save_emf = 1
  save_index = 1
  while (length(split_emfs) > 0) {
    curr_emf = names(split_emfs)[save_emf]
    curr_peaks = split_emfs[[save_emf]]
    for (iemf in seq_along(split_emfs)) {
      if (sum(curr_peaks %in% split_emfs[[iemf]]) == length(curr_peaks)) {
        curr_emf = unique(c(curr_emf, names(split_emfs)[iemf]))
        save_emf = unique(c(save_emf, iemf))
      }
    }
    emf_list[[save_index]] = curr_emf
    split_emfs = split_emfs[-save_emf]
    save_index = save_index + 1
    save_emf = 1
  }

  n_emf = purrr::map_int(emf_list, length)
  emf_list = emf_list[n_emf > 0]
  n_emf = n_emf[n_emf > 0]

  names(emf_list) = paste0("GEMF.", seq_along(emf_list))
  grouped_emf_mapping = data.frame(v1 = rep(names(emf_list), n_emf),
                                   v2 = unlist(emf_list, use.names = FALSE),
                                   stringsAsFactors = FALSE)
  names(grouped_emf_mapping) = c(grouping_var, emf_var)

  peak_info = dplyr::left_join(peak_info, grouped_emf_mapping, by = emf_var)
  peak_info
}

#' Count Carbon 13s
#'
#' Counts the number of carbons in the SMIRFE IMF formula's
#'
#' @param imf a character vector of isotopically resolved molecular formula's (IMFs)
#'
#' @export
#' @return numeric vector
get_c13_incorporation = function(imf){
  split_imfs = stringr::str_split(imf, ",")

  n_c13 = purrr::map_dbl(split_imfs, function(in_imf){
    #message(in_imf)
    if (sum(grepl("13C", in_imf)) == 1) {
      n_c = as.numeric(stringr::str_split_fixed(grep("13C", in_imf, value = TRUE), "C", 2))[2]
    } else {
      n_c = 0
    }
    n_c
  })

  n_c13
}
