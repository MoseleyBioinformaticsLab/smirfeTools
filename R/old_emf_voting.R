old_extract_emfs = function(chosen_emfs){
  all_samples = unique(unlist(purrr::map(chosen_emfs, ~ .x$Sample)))

  null_sample_matrix = matrix(as.character(NA), nrow = 1, ncol = length(all_samples))
  colnames(null_sample_matrix) = all_samples

  all_emfs = internal_map$map_function(chosen_emfs, function(peak_info){
  #all_emfs = purrr::map(seq(1, length(chosen_emfs)), function(i_emf){
    #message(i_emf)
    #peak_info = chosen_emfs[[i_emf]]

    split_emf = split(peak_info, peak_info$complete_EMF)

    emf_data = purrr::map(split_emf, function(in_emf){
      split_imf = split(in_emf, in_emf$complete_IMF)
      emf_matrix = purrr::map(split_imf, function(x){
        peak_matrix = null_sample_matrix
        peak_matrix[1, x$Sample] = x$Sample_Peak
        peak_matrix
      })

      emf_matrix = do.call(rbind, emf_matrix)

      emf_nap = purrr::map_df(split_imf, ~ dplyr::select(.x, complete_IMF, complete_EMF, NAP) %>% dplyr::slice(1))

      nap_order = order(emf_nap$NAP, decreasing = TRUE)

      list(peak_matrix = emf_matrix[nap_order, , drop = FALSE], peak_info = emf_nap[nap_order, , drop = FALSE])
    })

    peak_matrices = purrr::map(emf_data, ~ .x$peak_matrix)

    # there is an odd property of some of these EMFs. In some samples, one EMF will get assignments, but
    # the other EMF will not, leading to differences in the peak matrices. But if both EMFs have the same
    # number of isotopologues, then they should have both been seen in the same samples. This code
    # is an attempt to rectify that situation, by checking that either the same peaks were assigned the IMF
    # in some samples, or that no peaks were observed. If there is a sample with discrepant peaks observed
    # across the IMFs, then this will break and fall back to the original peak matrices
    n_row = purrr::map_int(peak_matrices, ~ nrow(.x))
    max_row = max(n_row)
    corresponded_matrices = peak_matrices

    for (irow in seq_len(max_row)) {
      use_matrices = irow <= n_row
      peak_compare = purrr::map_dfc(peak_matrices[use_matrices], ~ .x[irow, ]) %>% as.matrix()

      match_or_na = purrr::map_chr(seq_len(nrow(peak_compare)), function(in_peak){
        tmp_peak = peak_compare[in_peak, ]
        tmp_peak = tmp_peak[!is.na(tmp_peak)]
        unique_peak = length(unique(tmp_peak[!is.na(tmp_peak)])) == 1
        if (length(tmp_peak) == 0) {
          #warning("only NA!")
          return(as.character(NA))
        } else if (unique_peak) {
          return(tmp_peak[1])
        } else {
          #message("there was a mismatch!")
          return("FALSE")
        }
      })

      corresponded_matrices[use_matrices] = purrr::map(corresponded_matrices[use_matrices], function(.x){
        .x[irow, ] = match_or_na
        .x
      })

    }

    all_peaks = do.call(rbind, corresponded_matrices)

    if ("FALSE" %in% all_peaks) {
      #warning(paste0("there was a mismatch in ", i_emf))
      use_matrices = peak_matrices
    } else {
      use_matrices = corresponded_matrices
    }
    emf_data = purrr::map(seq_len(length(emf_data)), function(idata){
      emf_data[[idata]]$corresponded_matrix = use_matrices[[idata]]
      emf_data[[idata]]
    })

    names(emf_data) = names(split_emf)

    emf_data

  })

  names(all_emfs) = names(chosen_emfs)
  semf_emf = purrr::map2_df(all_emfs, names(all_emfs), function(.x, .y){
    data.frame(semf = .y, emf = names(.x), stringsAsFactors = FALSE)
  })

  if (sum(duplicated(semf_emf$emf)) > 0) {
    warning("Duplicate sudo EMF to complete EMF mappings detected!")
  }

  all_emfs
}

#' add peak location and intensity
#'
#' Given a set of EMF peak matrices and a data.frame of various locations and
#' intensities, associates a matrix of each of location and intensity directly
#' to the EMF peak matrix. Can be used to extract alternatives to the default
#' of Height and ObservedMZ.
#'
#' @param emfs list of EMFs
#' @param peak_data data.frame of peak data
#' @param location which variable to use for location (unquoted)
#' @param intensity which variable to use for intensity (unquoted)
#'
#' @export
#' @return list
old_add_location_intensity = function(emfs, peak_data, location = ObservedMZ,
                                  intensity = Height){

  peak_data = dplyr::mutate(peak_data, location = {{location}}, intensity = {{intensity}})
  peak_location = dplyr::select(peak_data, Sample_Peak, location)
  peak_intensity = dplyr::select(peak_data, Sample_Peak, intensity)

  peak_intensity = rbind(peak_intensity, data.frame(Sample_Peak = "0", intensity = NA, stringsAsFactors = FALSE))

  peak_intensity_matrix = matrix(peak_intensity$intensity, nrow = nrow(peak_intensity), ncol = 1)
  rownames(peak_intensity_matrix) = peak_intensity$Sample_Peak

  peak_location = rbind(peak_location, data.frame(Sample_Peak = "0", location = NA,
                                      stringsAsFactors = FALSE))
  peak_location_matrix = matrix(peak_location$location, nrow = nrow(peak_location), ncol = 1)
  rownames(peak_location_matrix) = peak_location$Sample_Peak

  all_samples = colnames(emfs[[1]][[1]]$peak_matrix)
  numeric_matrix = matrix(NA, nrow = 1, ncol = length(all_samples))
  colnames(numeric_matrix) = all_samples
  emfs_2 = internal_map$map_function(emfs, function(in_emf){
    purrr::map(in_emf, function(group_emf){
      tmp_intensity = purrr::map(seq(1, nrow(group_emf$corresponded_matrix)), function(in_row){
        tmp_data = numeric_matrix
        tmp_peaks = group_emf$corresponded_matrix[in_row, ]
        tmp_peaks[is.na(tmp_peaks)] = "0"
        tmp_data[1, ] = peak_intensity_matrix[tmp_peaks, 1]

      })
      out_intensity = do.call(rbind, tmp_intensity)
      colnames(out_intensity) = all_samples
      group_emf$intensity = out_intensity

      tmp_location = purrr::map(seq(1, nrow(group_emf$corresponded_matrix)), function(in_row){
        tmp_data = numeric_matrix
        tmp_peaks = group_emf$corresponded_matrix[in_row, ]
        tmp_peaks[is.na(tmp_peaks)] = "0"
        tmp_data[1, ] = peak_location_matrix[tmp_peaks, 1]

      })
      out_location = do.call(rbind, tmp_location)
      colnames(out_location) = all_samples
      group_emf$location = out_location

      group_emf
    })
  })

}

#' extract IMF or EMF level data
#'
#' This function is deprecated, and kept around merely for use in some old code.
#' Extract either IMF or EMF level data from the extracted EMFs. Which is controlled by the `by` argument.
#'
#' @param emfs list of EMFs
#' @param by `EMF` (default) or `IMF`
#'
#' @details For `EMF`, will be a list corresponding to each EMF from the sudo EMFs, for `IMF`, will be
#'  matrices and a data.frame.
#'
#' @return list with `intensity`, `location`, and `info`
#' @export
#'
old_extract_imf_emf_data = function(emfs, by = "EMF"){
  .Deprecated("extract_imf_emf_data")
  all_compared = internal_map$map_function(emfs, old_compare_extract_emfs, by = by)

  if (is.null(names(all_compared))) {
    names(all_compared) = paste0("SEMF_", seq_along(all_compared))
  }

  if ("EMF" %in% by) {
    all_intensity = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_intensity = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, ".", emf_id, ".", names(in_emf$info))
        out_intensity = in_emf$intensity
        rownames(out_intensity) = peak_id
        out_intensity
      })
      names(semf_intensity) = paste0(semf_id, ".", names(in_semf))
      semf_intensity

    }) %>% purrr::flatten()

    all_location = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_location = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, ".", emf_id, ".", names(in_emf$info))
        out_location = in_emf$location
        rownames(out_location) = peak_id
        out_location
      })
      names(semf_location) = paste0(semf_id, ".", names(in_semf))
      semf_location

    }) %>% purrr::flatten()

    all_info = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_info = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, ".", emf_id, ".", names(in_emf$info))
        info_df = purrr::map2_dfr(in_emf$info, peak_id, function(in_info, pid){
          in_info$PeakID = pid
          in_info$sudo_EMF = semf_id
          in_info$sudo_group_EMF = paste0(emf_id)
          in_info
        })
      })
      names(semf_info) = paste0(semf_id, ".", names(in_semf))
      semf_info
    }) %>% purrr::flatten()

    all_peaks = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_peaks = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, ".", emf_id, ".", names(in_emf$info))
        out_peaks = in_emf$peaks
        rownames(out_peaks) = peak_id
        out_peaks
      })
      names(semf_peaks) = paste0(semf_id, ".", names(in_semf))
      semf_peaks

    }) %>% purrr::flatten()

  } else {
    all_intensity = purrr::map2(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,".", names(.x$info))
      out_intensity = .x$intensity
      n_intensity = nrow(out_intensity)
      #message(paste0(.y, " ", n_height))
      rownames(out_intensity) = peak_id
      out_intensity
    }) %>% do.call(rbind, .)
    all_location = purrr::map2(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,".", names(.x$info))
      out_location = .x$location
      rownames(out_location) = peak_id
      out_location
    }) %>% do.call(rbind, .)
    all_info = purrr::map2_dfr(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,".", names(.x$info))
      info = .x$info
      n_info = purrr::map_int(info, nrow)
      info_df = purrr::map2_dfr(info, peak_id, function(tmp_info, tmp_peak){
        tmp_info$PeakID = tmp_peak
        tmp_info$sudo_EMF = .y
        tmp_info
      })
      info_df
    })
  }
  return(list(intensity = all_intensity,
              location = all_location,
              peaks = all_peaks,
              info = all_info))
}

old_compare_extract_emfs = function(emf, by = "EMF"){

  get_imfs = function(emf_group){
    split_groups = split(emf_group, emf_group$group_imfs)
    names(split_groups) = paste0("IMF_", seq_along(split_groups))
    single_index = purrr::map_df(split_groups, ~ .x[1, ])

    out_intensity = compare_intensity[single_index$row_index, , drop = FALSE]
    out_location = compare_location[single_index$row_index, , drop = FALSE]
    out_peaks = compare_peaks[single_index$row_index, , drop = FALSE]

    out_info = purrr::map(split_groups, ~ compare_info[.x$row_index, ])


    list(intensity = out_intensity, location = out_location, info = out_info,
         peaks = out_peaks)
  }

  compare_intensity = purrr::map(emf, ~ .x$intensity) %>% do.call(rbind, .)
  compare_info = purrr::map_df(emf, ~ .x$peak_info)
  compare_location = purrr::map(emf, ~ .x$location) %>% do.call(rbind, .)
  compare_peaks = purrr::map(emf, ~ .x$corresponded_matrix) %>% do.call(rbind, .)

  group_imfs = vector("integer", length = nrow(compare_intensity))

  i_group = 1
  while (sum(group_imfs == 0) > 0) {
    zero_locs = which(group_imfs == 0)
    master_loc = which(group_imfs == 0)[1]

    for (iloc in zero_locs) {
      if (isTRUE(all.equal(compare_intensity[master_loc, ], compare_intensity[iloc, ]))) {
        group_imfs[iloc] = i_group
      }
    }
    i_group = i_group + 1
  }

  df_groups = data.frame(group_imfs = group_imfs, row_index = seq(1, length(group_imfs)))
  df_groups$complete_EMF = compare_info$complete_EMF[df_groups$row_index]

  if ("EMF" %in% by) {
    split_groups = split(df_groups, df_groups$complete_EMF)
    group_emfs = vector("integer", length = length(split_groups))

    i_group = 1
    while (sum(group_emfs == 0) > 0) {
      zero_locs = which(group_emfs == 0)
      master_loc = which(group_emfs == 0)[1]

      for (iloc in zero_locs) {
        if (isTRUE(all.equal(split_groups[[master_loc]]$group_imfs, split_groups[[iloc]]$group_imfs))) {
          group_emfs[iloc] = i_group
        }
      }
      i_group = i_group + 1
    }

    df_groups2 = purrr::map2_dfr(split_groups, group_emfs, function(.x, .y){
      .x$group_emfs = .y
      .x
    })
    split_emfs = split(df_groups2, df_groups2$group_emfs)
    names(split_emfs) = paste0("EMF_", seq_along(split_emfs))

    out_data = purrr::map(split_emfs, get_imfs)

  } else {
    out_data = get_imfs(df_groups)
  }

  out_data

}

#' extract EMF propensity
#'
#' If one wants to be able to filter EMF based data for those EMFs that show up in a particular number
#' of samples, one needs to be able to count which samples the EMF appears in. This function extracts
#' creates a matrix where ones note the samples that contained the EMF.
#'
#' @param emf_intensity list of EMF intensities
#'
#' @export
#' @return matrix
#'
old_extract_emf_propensity = function(emf_intensity){
  presence_matrix = purrr::map(emf_intensity, ~ as.integer(colSums(.x) > 0)) %>% do.call(rbind, .)
  rownames(presence_matrix) = names(emf_intensity)
  colnames(presence_matrix) = colnames(emf_intensity[[1]])
  presence_matrix
}
