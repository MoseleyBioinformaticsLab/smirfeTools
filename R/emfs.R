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


#' Decide EMF Peaks
#'
#' A given grouped EMF might have several possible peaks and complete EMF possibilities.
#' This function attempts to narrow those down by considering if there is
#' a lipid category prediction, and then the difference between the maximum
#' number of C13 incorporations observed and the actual number of peaks
#' in the data for that EMF. Only those EMFs that have the same number of
#' peaks are returned.
#'
#' @param emf_peaks *all* peaks for a grouped EMF, as a data.frame. See Details.
#'
#' @return data.frame
#' @export
decide_emf_peaks = function(emf_peaks){
  split_by_emf = split(emf_peaks, emf_peaks$complete_EMF)

  split_by_emf = purrr::map(split_by_emf, function(x){
    x = unique(dplyr::select(x, -isotopologue_EMF, -SortedFormula))
    x = dplyr::filter(x, !(is.na(PredictedCategories) & is.na(Categories)))
    x
  })

  n_res = purrr::map_int(split_by_emf, nrow)
  n_class = purrr::map_int(split_by_emf, ~ length(unique(.x$PredictedCategories)))

  split_by_emf = split_by_emf[(n_res > 0) | (n_class == 1)]

  if (length(split_by_emf) == 0) {
    return(emf_peaks[0, ])
  }

  has_k = stringr::str_detect(names(split_by_emf), "K")

  check_n_k_isotopes = function(in_emf){
    in_emf = dplyr::mutate(in_emf, k_iso = stringr::str_match(complete_IMF, "[[:digit:]]{2}K")[,1])

    split_c13 = split(in_emf, in_emf$C13)
    count_k = purrr::map_int(split_c13, nrow)
    if (any(count_k > 1)) {
      out_emf = split(in_emf, in_emf$k_iso)
      out_emf = purrr::map(out_emf, function(x){
        x$complete_EMF = paste0(x$complete_EMF, ".", x$k_iso)
        x$k_iso = NULL
        x
      })
    } else {
      in_emf$k_iso = NULL
      out_emf = list(in_emf)
      #names(out_emf) = in_emf[1, "complete_EMF"]
    }

    out_emf
  }

  if (sum(has_k) >= 1) {
    tmp_k = purrr::map(split_by_emf[has_k], check_n_k_isotopes)
    tmp_k = unlist(tmp_k, recursive = FALSE)
    split_by_emf = c(split_by_emf[!has_k], tmp_k)
  }

  if (length(split_by_emf) > 1) {
    emf_imf_info = purrr::imap_dfr(split_by_emf, function(.x, .y){
      max_imf = max(.x$C13) + 1
      n_imf = length(unique(.x$peak))

      dist_matrix = matrix(0L, nrow = 2, ncol = max_imf)
      dist_matrix[1, ] = 1
      dist_matrix[2, .x$C13+1] = 1
      diff_imf = dist(dist_matrix)
      data.frame(n = n_imf, diff = diff_imf[1], emf = .y, stringsAsFactors = FALSE)
    })

    emf_imf_info = dplyr::filter(emf_imf_info, diff == min(diff))

    use_emfs = dplyr::filter(emf_imf_info, n == max(n))

    keep_emfs = purrr::map_df(split_by_emf[use_emfs$emf], ~ .x)

    n_peak = length(unique(keep_emfs$peak))

    if (n_peak == use_emfs[1, "n"]) {
      out_emfs = keep_emfs
    } else {
      out_emfs = keep_emfs[0, ]
    }
  } else {
    out_emfs = split_by_emf
  }

  out_emfs
}

#' Cleanup NA Categories
#'
#' After choosing a set of peaks for a grouped EMF, it may be that there are
#' predicted categories that have an NA match in the Category based on formula
#' matching. In this particular case, we actually want to keep just the ones
#' that are NOT NA in both. This function does that. It also checks that we
#' have the same number of peaks as went in. Otherwise, it just returns the
#' original data.
#'
#' @param grouped_emf a data.frame of peaks for the grouped emf
#'
#' @export
#' @return data.frame
cleanup_na_categories = function(grouped_emf){
  pred_na = dplyr::filter(grouped_emf, !is.na(PredictedCategories), !is.na(Categories))

  if (length(unique(pred_na$peak)) == length(unique(grouped_emf$peak))) {
    return(pred_na)
  } else {
    return(grouped_emf)
  }

}