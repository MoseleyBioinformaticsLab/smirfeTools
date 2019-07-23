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
group_emfs_by_peaks = function(peak_info, peak_var = "peak", emf_var = "complete_EMF",
                               grouping_var = "grouped_EMF"){
  split_emfs = split(peak_info[, peak_var], peak_info[, emf_var])

  emf_list = vector("list", length(split_emfs))

  # this finds the formulas that share identical sets of peaks
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

  peak_info_gemf = dplyr::left_join(peak_info, grouped_emf_mapping, by = emf_var)

  # now we need to find those grouped EMFs that *share* peaks. This is a no-no, because
  # peaks can't be part of two EMFs. The solution right now is to keep the EMF with the
  # most peaks, and remove the other EMFs.
  split_gemfs = split(peak_info_gemf, peak_info_gemf[, grouping_var])
  gemfs_df = purrr::map_df(split_gemfs, ~ unique(.x[, c(peak_var, grouping_var)]))
  names(gemfs_df) = c("peak", "group")

  gemfs_merge = dplyr::left_join(gemfs_df, gemfs_df, by = "peak", suffix = c(".1", ".2"))

  gemfs_notequal = gemfs_merge[(gemfs_merge$group.1 != gemfs_merge$group.2), ]

  if (nrow(gemfs_notequal) > 0) {

    gemfs_compare = vector("list", nrow(gemfs_notequal))
    compare_loc = 1
    while (nrow(gemfs_notequal) > 0) {
      init_gemfs = unique(unlist(gemfs_notequal[1, c("group.1", "group.2")], use.names = FALSE))
      more_gemfs = unique(unlist(dplyr::filter(gemfs_notequal, (group.1 %in% init_gemfs) | (group.2 %in% init_gemfs))[, c("group.1", "group.2")], use.names = FALSE))
      while (!all(more_gemfs %in% init_gemfs)) {
        init_gemfs = more_gemfs
        more_gemfs = unique(unlist(dplyr::filter(gemfs_notequal, (group.1 %in% init_gemfs) | (group.2 %in% init_gemfs))[, c("group.1", "group.2")], use.names = FALSE))
      }
      gemfs_compare[[compare_loc]] = more_gemfs
      compare_loc = compare_loc + 1
      gemfs_notequal = gemfs_notequal[!((gemfs_notequal$group.1 %in% more_gemfs) | (gemfs_notequal$group.2 %in% more_gemfs)), ]
    }

    n_merge = purrr::map_int(gemfs_compare, length)
    gemfs_compare = gemfs_compare[n_merge > 0]

    count_peaks_lost = purrr::map_int(gemfs_compare, function(in_compare){
      tmp_matrix = create_peaks_2_emf(split_gemfs[in_compare])
      calculate_peaks_lost(tmp_matrix)
    })

    gemfs_df_npeak = purrr::map_int(split(gemfs_df$peak, gemfs_df$group), ~ length(.x))

    gemfs_compare_keep = purrr::map_chr(gemfs_compare, function(compare_group){
      gemfs_count = gemfs_df_npeak[compare_group]
      names(gemfs_count)[which.max(gemfs_count)]
    })

    not_equal = unique(unlist(gemfs_compare, use.names = FALSE))

    all_emfs = names(emf_list)
    keep_emfs = unique(c(setdiff(all_emfs, not_equal), gemfs_compare_keep))
    emf_list2 = emf_list[keep_emfs]
    names(emf_list2) = paste0("GEMF.", seq_along(emf_list2))

    n_emf2 = purrr::map_int(emf_list2, length)

    grouped_emf_mapping2 = data.frame(v1 = rep(names(emf_list2), n_emf2),
                                     v2 = unlist(emf_list2, use.names = FALSE),
                                     stringsAsFactors = FALSE)
    names(grouped_emf_mapping2) = c(grouping_var, emf_var)

    peak_info_gemf2 = dplyr::left_join(peak_info, grouped_emf_mapping2, by = emf_var)
    peak_info_gemf2 = peak_info_gemf2[!is.na(peak_info_gemf2[, grouping_var]), ]

    peak_info_gemf = peak_info_gemf2

  }
  peak_info_gemf
}

create_peaks_2_emf = function(emf_list){
  emf_peaks = purrr::map(emf_list, ~ unique(.x$peak))
  all_peaks = unique(unlist(emf_peaks, use.names = FALSE))

  peak_matrix = matrix(0, nrow = length(all_peaks), ncol = length(emf_list))
  colnames(peak_matrix) = names(emf_list)
  rownames(peak_matrix) = all_peaks

  for (i_emf in names(emf_list)) {
    peak_matrix[unique(emf_list[[i_emf]]$peak), i_emf] = 1
  }
  peak_matrix
}

calculate_peaks_lost = function(peak_matrix){
  biggest_emf = which.max(colSums(peak_matrix))

  missing_peaks = (rowSums(peak_matrix) > 0) & (peak_matrix[, biggest_emf] == 0)
  n_peak = sum(missing_peaks)
  n_peak
}

find_non_zeros = function(peak_matrix){
  purrr::map_dbl(seq(1, nrow(peak_matrix)), function(x){
    tmp_loc = which(peak_matrix[x, ] > 0)
    if (length(tmp_loc) == 0) {
      tmp_loc = 0
    }
    tmp_loc
  })
}

permuate_peak_matrix = function(peak_matrix){
  peak_sums = rowSums(peak_matrix)
  peak_matrix = peak_matrix[order(peak_sums), ]
  peak_sums = peak_sums[order(peak_sums)]

  modify_rows = which(peak_sums > 1)

  mod_locs = purrr::map(modify_rows, ~ which(mod_matrix[.x, ] > 0))

  mod_matrix = peak_matrix

  # set all the multi EMF peaks to 0 for a baseline
  for (irow in names(mod_locs)) {
    mod_matrix[irow, ] = 0
  }

  score_zero = score_function(mod_matrix)
  zero_permutation = find_non_zeros(mod_matrix)

  # introduce each row, keeping the others zero
  zero_matrix = mod_matrix
  score_permutations = list()
  for (irow in names(mod_locs)) {
    use_locs = mod_locs[[irow]]
    for (icol in use_locs) {
      zero_matrix[irow, icol] = 1
      score_permutations = c(score_permuations, list(score = score_function(zero_matrix),
                                                     permutation = find_non_zeros(zero_matrix)))
      zero_matrix[irow, ] = 0
    }
  }
  # after this, we need to hold one row constant, and then permute all the other rows above and below it

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
    out_emfs = split_by_emf[[1]]
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
  pred_na = dplyr::filter(grouped_emf, !is.na(PredictedCategories))

  if ((length(unique(pred_na$peak)) == length(unique(grouped_emf$peak))) || (nrow(pred_na) == 0)) {
    return(pred_na)
  } else {
    return(grouped_emf)
  }

}

#' keep EMFs with matching features
#'
#' Given a set of features of interest, and a set of decided EMFs, goes through
#' and returns the full EMFs that contain any of the features of interest. Also
#' attempts to return as few classes of EMFs as possible.
#'
#' @param features a set of interesting features
#' @param emf_info data.frame of EMFs, including the lipid class information
#'
#' @export
#' @return list
keep_emfs_with_features = function(features, emf_info){
  peak_info = dplyr::filter(emf_info, peak %in% features)
  peak_gemfs = dplyr::filter(emf_info, grouped_EMF %in% unique(peak_info$grouped_EMF))

  peak_gemfs_split = base::split(peak_gemfs, peak_gemfs$complete_EMF)

  peak_gemfs_split = purrr::map(peak_gemfs_split, cleanup_na_categories)

  peak_gemfs_ncategory = purrr::map_int(peak_gemfs_split, ~ length(unique(.x$PredictedCategories)))

  peak_gemfs_split = peak_gemfs_split[peak_gemfs_ncategory == 1]
}

#' Get EMF sums
#'
#' @param split_emfs list of EMF data.frames
#' @param feature_intensities matrix of feature (rows) by sample (column) intensities
#' @param feature_mz matrix of M/Z values
#' @param multiply_C13 should the intensities be multiplied by the C13 incorporation? (FALSE)
#'
#' @importFrom dplyr "%>%"
#' @export
#' @return list
get_emf_intensity_sums = function(split_emfs, feature_intensities, feature_mz, multiply_C13 = FALSE){

  emf_data = list(sums = purrr::map(split_emfs, function(in_emf){
    if (multiply_C13){
      c13_vals = in_emf$C13
    } else {
      c13_vals = rep(1, nrow(in_emf))
    }
    colSums(feature_intensities[in_emf$peak, , drop = FALSE] * c13_vals)
  }) %>% do.call(base::rbind, .),
  classes = purrr::map_df(split_emfs, function(in_emf){
    dplyr::select(in_emf, complete_EMF, PredictedClasses,
                  PredictedCategories, Categories, Classes) %>%
      dplyr::sample_n(., 1)}) %>%
    dplyr::mutate(., mz = purrr::map_dbl(split_emfs, function(in_emf){
                    use_peak = in_emf$peak[in_emf$C13 == min(in_emf$C13)]
                    min(feature_mz[use_peak, ], na.rm = TRUE)
                  }))
  )
  emf_data
}

#' Extract Single M/Z
#'
#' Given a list of EMF M/Z, we would like to get a single M/Z for each EMF,
#' as by definition an EMF often has at least two IMFs.
#'
#' @param emf_list_mz the list of EMF M/Zs
#'
#' @export
#'
extract_emf_mz = function(emf_list_mz){
  single_mz = function(matrix_mz, emf_id){
    mean_samples = rowMeans(matrix_mz, na.rm = TRUE)
    data.frame(emf = emf_id, mz = mean_samples[1], stringsAsFactors = FALSE)
  }

  emf_mz = purrr::imap_dfr(emf_list_mz, single_mz)
  emf_mz
}
