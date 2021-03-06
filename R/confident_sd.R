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
#' @return data.frame with sd values
find_confident_frequency_sd = function(assigned_data,
                             low_evalue_cutoff = 0.1,
                             low_mz_cutoff = 600,
                             remove_elements = "S"){

  scan_level_frequency = purrr::map(assigned_data, ~ .x$scan_level$ObservedFrequency)
  scan_level_names = unlist(purrr::map(scan_level_frequency, ~ names(.x)))
  scan_level_frequency = unlist(scan_level_frequency, recursive = FALSE, use.names = FALSE)
  names(scan_level_frequency) = scan_level_names

  sample_peak = purrr::map_df(assigned_data, ~ unique(.x$data[, c("Sample", "Sample_Peak")]))

  confident_emfs = purrr::map(assigned_data, function(in_assign){
    assignments = in_assign$assignments
    #message(paste0(in_assign, "  ", .x$sample))
    low_e_mz = dplyr::filter(assignments, (e_value <= low_evalue_cutoff) &
                               (ObservedMZ <= low_mz_cutoff) & (!grepl("S", complete_EMF)))

    if (length(unique(low_e_mz$Sample_Peak)) >= 20) {
      return(get_sample_emfs(low_e_mz, in_assign$sample, evalue_cutoff = low_evalue_cutoff, use_corroborating = FALSE))
    } else {
      return(NULL)
    }

  })

  confident_emfs = confident_emfs[!purrr::map_lgl(confident_emfs, is.null)]
  confident_gemf_emf_mapping = internal_map$map_function(confident_emfs, function(x){
    purrr::map_df(x, ~ unique(dplyr::select(.x, grouped_EMF, complete_EMF)))
  })
  confident_gemf_emf_mapping = do.call(rbind, confident_gemf_emf_mapping)

  confident_sudo_emfs = create_sudo_emfs(confident_gemf_emf_mapping)

  n_emf_confident = purrr::map_int(confident_sudo_emfs, ~length(unique(.x$grouped_EMF)))
  confident_sudo_emfs = confident_sudo_emfs[n_emf_confident > 1]
  confident_all_gemfs = unlist(confident_emfs, recursive = FALSE, use.names = FALSE)
  names(confident_all_gemfs) = purrr::map_chr(confident_all_gemfs, ~ .x$grouped_EMF[1])

  # sd_information = internal_map$map_function(confident_sudo_emfs, function(in_sudo){
  #   calculate_confident_sd(confident_all_gemfs[unique(in_sudo$grouped_EMF)], scan_level_frequency, sample_peak)
  # })

  sd_information = internal_map$map_function(names(confident_sudo_emfs), function(sudo_id){
    #message(sudo_id)
    "!DEBUG `sudo_id`"
    in_sudo = confident_sudo_emfs[[sudo_id]]
    calculate_confident_sd(confident_all_gemfs[unique(in_sudo$grouped_EMF)], scan_level_frequency, sample_peak)
  })
  names(sd_information) = names(confident_sudo_emfs)

  sd_df = purrr::imap_dfr(sd_information, function(.x, .y){
    data.frame(semf = .y, sd = .x, stringsAsFactors = FALSE)
  })

  sd_df

}

calculate_confident_sd = function(emf_data, scan_level_frequency, sample_peak){
  peaks_2_emf_imf = purrr::map_df(emf_data,  ~ unique(.x[, c("complete_EMF", "complete_IMF", "Sample_Peak")]))

  split_peaks = split(peaks_2_emf_imf$Sample_Peak, peaks_2_emf_imf$complete_IMF)

  split_peaks = purrr::map(split_peaks, unique)
  split_peaks = unique(split_peaks)

  n_peak = purrr::map_int(split_peaks, length)
  split_peaks = split_peaks[order(n_peak, decreasing = TRUE)]


  # this code is trying to do the same thing as for reducing duplicates in
  # sudo emfs, basically testing if there is more than 50% overlap of peaks
  # between IMFs, and if so, merges them.
  if (length(split_peaks) > 1) {
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

  } else {
    merged_peaks = split_peaks
  }


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


#' create M/Z differences
#'
#' Given a list of coefficient data from multiple samples and a previously
#' derived `frequency_offset`, this generates a data.frame of M/Z locations
#' and differences.
#'
#' @param coefficient_info list of coefficients from `get_sample_coefficients`
#' @param frequency_offset the offset to use, based on `find_confident_frequency_sd`
#' @param mz_values a vector of M/Z values to create offsets for
#'
#' @export
#' @return data.frame with Index and Value
create_mz_diffs = function(coefficient_info, frequency_offset = NULL, mz_values = seq(150, 1600, 0.5)){
  if (is.null(frequency_offset)) {
    stop("frequency_offset is NULL, please supply a value!")
  }

  freq_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$frequency_coefficients
  }) %>% do.call(rbind, .) %>% apply(., 2, median)

  mz_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$mz_coefficients
  }) %>% do.call(rbind, .) %>% apply(., 2, median)

  freq_description = coefficient_info[[1]]$coefficients$frequency_fit_description
  mz_description = coefficient_info[[1]]$coefficients$mz_fit_description

  mz_2_freq = predict_exponentials(mz_values, freq_coefficients, freq_description)

  frequency_diff = mz_2_freq + frequency_offset

  mz_init = predict_exponentials(mz_2_freq, mz_coefficients, mz_description)
  mz_shift = predict_exponentials(frequency_diff, mz_coefficients, mz_description)
  mz_diff = abs(mz_shift - mz_init)

  data.frame(Index = mz_values, Value = mz_diff)

}

#' get coefficient data
#'
#' Given a set of zip files, find the associated json metadata files, and extract the
#' coefficient information from them.
#'
#' @param zip_file_list the list of zip files
#'
#' @export
#' @return list of lists
extract_coefficient_data = function(zip_file_list){
  json_files = gsub(".zip$", ".json", zip_file_list)

  exist_json = file.exists(json_files)

  if (sum(!exist_json) > 0) {
    warning_message = paste0(sum(!exist_json), " json files of ", length(json_files),
                             " couldn't be found!")
    warning(warning_message)
  }

  internal_map$map_function(json_files[exist_json], function(use_file){
    list(sample = gsub(".json$", "", basename(use_file)),
         coefficients = jsonlite::fromJSON(use_file)$peak$frequency_mz)
  })
}
