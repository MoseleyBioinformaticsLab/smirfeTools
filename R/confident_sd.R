#' find confident frequency sd
#'
#' Given a list of EMF assignment data from SMIRFE, attempts to find *confident*
#' EMF assignments, and use those combined with scan level frequency information
#' to determine a suitable cutoff across sample assignments.
#'
#' @param assigned_data list of assignment data
#' @param low_evalue_cutoff upper bound for what is considered a *confident* assignment
#' @param low_mz_cutoff the upper bound for *confident* assignments in M/Z
#'
#' @export
#' @return
find_confident_frequency_sd = function(assigned_data,
                             low_evalue_cutoff = 0.1,
                             low_mz_cutoff = 600,
                             remove_elements = "S"){

  scan_level_frequency = purrr::map(assigned_data, ~ .x$scan_level$ObservedFrequency)
  scan_level_names = unlist(purrr::map(scan_level_frequency, ~ names(.x)))
  scan_level_frequency = unlist(scan_level_frequency, recursive = FALSE, use.names = FALSE)
  names(scan_level_frequency) = scan_level_names

  sample_peak = purrr::map_df(assigned_data, ~ unique(.x$data[, c("Sample", "Sample_Peak")]))

  confident_emfs = internal_map$map_function(seq(1, length(assigned_data)), function(in_assign){
    .x = assigned_data[[in_assign]]
    #message(paste0(in_assign, "  ", .x$sample))
    peak_mz = dplyr::filter(.x$data, Measurement %in% "ObservedMZ")
    low_mz_peaks = dplyr::filter(peak_mz, Value <= low_mz_cutoff) %>% dplyr::pull(Sample_Peak)

    tmp_assign = dplyr::filter(.x$assignments, !grepl(remove_elements, complete_EMF))
    assign_evalues = dplyr::filter(tmp_assign, Type %in% "e_value")
    low_e_peaks = dplyr::filter(assign_evalues, as.numeric(Assignment_Data) <= low_evalue_cutoff) %>% dplyr::pull(Sample_Peak)

    if (length(base::intersect(low_mz_peaks, low_e_peaks)) >= 20) {
      return(get_sample_emfs(tmp_assign, .x$sample, evalue_cutoff = low_evalue_cutoff))
    } else {
      return(NULL)
    }

  })

  confident_emfs = confident_emfs[!purrr::map_lgl(confident_emfs, is.null)]
  confident_gemf_emf_mapping = internal_map$map_function(confident_emfs, function(x){
    purrr::map2_dfr(x$grouped_emf, names(x$grouped_emf), function(.x, .y){
      data.frame(grouped_emf = .y, complete_EMF = .x$complete_EMF, stringsAsFactors = FALSE)
    })
  })
  confident_gemf_emf_mapping = do.call(rbind, confident_gemf_emf_mapping)

  confident_sudo_emfs = create_sudo_emfs(confident_gemf_emf_mapping)

  n_emf_confident = purrr::map_int(confident_sudo_emfs, ~length(unique(.x$grouped_emf)))
  confident_sudo_emfs = confident_sudo_emfs[n_emf_confident > 1]
  confident_all_gemfs = unlist(purrr::map(confident_emfs, "grouped_emf"), recursive = FALSE)

  sd_information = internal_map$map_function(confident_sudo_emfs, function(in_sudo){
    calculate_confident_sd(confident_all_gemfs[unique(in_sudo$grouped_emf)], scan_level_frequency, sample_peak)
  })

  sd_df = purrr::imap_dfr(sd_information, function(.x, .y){
    data.frame(semf = .y, sd = .x, stringsAsFactors = FALSE)
  })

  sd_df

}

calculate_confident_sd = function(emf_data, scan_level_frequency, sample_peak){
  peaks_2_emf_imf = purrr::map_df(emf_data, function(in_data){
    purrr::map_df(in_data$peak_info, ~ unique(.x[, c("complete_EMF", "complete_IMF", "Sample_Peak")]))
  })

  split_peaks = split(peaks_2_emf_imf$Sample_Peak, peaks_2_emf_imf$complete_IMF)

  split_peaks = purrr::map(split_peaks, unique)
  split_peaks = unique(split_peaks)

  n_peak = purrr::map_int(split_peaks, length)
  split_peaks = split_peaks[order(n_peak, decreasing = TRUE)]

  merged_peaks = vector("list", length(split_peaks))
  combined_data = merged_peaks
  is_grabbable = rep(TRUE, length(split_peaks))

  for (isplit in seq(1, length(split_peaks) - 1)) {
    if (is_grabbable[isplit]) {
      merged_peaks[[isplit]] = split_peaks[[isplit]]
      is_grabbable[isplit] = FALSE
    }

    for (jsplit in seq(isplit+1, length(split_peaks))) {
      #message(jsplit)
      if (is_grabbable[jsplit]) {
        inter_ratio = length(intersect(split_peaks[[jsplit]], merged_peaks[[isplit]])) / length(split_peaks[[jsplit]])
        if (inter_ratio >= 0.5) {
          merged_peaks[[isplit]] = union(split_peaks[[jsplit]], merged_peaks[[isplit]])
          is_grabbable[jsplit] = FALSE
          combined_data[[isplit]] = unique(c(combined_data[[isplit]], isplit, jsplit))
        }
      }
    }
  }

  merged_peaks = merged_peaks[!(purrr::map_lgl(merged_peaks, is.null))]

  merged_peaks = purrr::map(merged_peaks, function(in_peaks){
    tmp_sample = dplyr::filter(sample_peak, Sample_Peak %in% in_peaks)
    if (nrow(tmp_sample) != length(in_peaks)) {
      dup_samples = tmp_sample$Sample[duplicated(tmp_sample$Sample)]
      tmp_sample = dplyr::filter(tmp_sample, Sample %in% dup_samples)
      message("found some!")
      return(tmp_sample$Sample_Peak)
    } else {
      return(in_peaks)
    }
  })

  sd_data = purrr::map_dbl(merged_peaks, function(in_peaks){
    in_frequency = unlist(scan_level_frequency[in_peaks])
    sd(in_frequency, na.rm = TRUE)
  })
  sd_data
}