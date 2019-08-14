#' extract sample level elemental molecular formulas
#'
#' given the assignment data.frame, group and extract EMF level information
#'
#' @param sample_assignments a sample level assignment data.frame from `read_smirfe_assignment`
#' @param sample_id which sample is it
#' @param evalue_cutoff what value should be used to exclude EMFs? (default is 0.98)
#' @param use_corroborating should other adducts be added for corroborating evidence (default is `TRUE`)
#'
#' @return list
#'
#' @export
get_sample_emfs = function(sample_assignments, sample_id, evalue_cutoff = 0.98, use_corroborating = TRUE){
  sample_assignments = dplyr::filter(sample_assignments, !grepl("S", complete_EMF))
  # add the adduct information into the EMF information
  sample_assignments = dplyr::mutate(sample_assignments, emf_Adduct = paste0(complete_EMF, ".", adduct_IMF))
  e_values = dplyr::filter(sample_assignments, Type %in% "e_value") %>% dplyr::mutate(e_value = as.numeric(Assignment_Data), imf_peak = paste0(complete_IMF, "_", adduct_IMF, "_", PeakID))
  clique_size = dplyr::filter(sample_assignments, Type %in% "clique_size") %>% dplyr::mutate(clique_size = as.integer(Assignment_Data), imf_peak = paste0(complete_IMF, "_", adduct_IMF, "_", PeakID))
  evalue_clique_size = dplyr::left_join(e_values, clique_size[, c("imf_peak", "clique_size")], by = "imf_peak")


  # split the peaks by which EMF they were assigned to, and then for each
  # EMF spit out the set of peaks that were assigned that EMF, both in an
  # easy to compare character form, and the actual list
  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$complete_EMF) %>%
    purrr::map2_dfr(., names(.), function(.x, .y){
      peaks = unique(.x)
      tmp_frame = data.frame(complete_EMF = .y, PeakID_chr = paste(peaks, collapse = ","), stringsAsFactors = FALSE)
      tmp_frame$PeakID = list(peaks)
      tmp_frame
    })

  # now we can easily find those EMFs that share the same set of peaks to create
  # `grouped_EMFs`
  grouped_emf = split(peaks_by_emf, peaks_by_emf$PeakID_chr)

  # for each grouped_EMF, get out the peak info, the e-value
  # that corresponds to the number of peaks in the clique,
  # filter to those EMF assignments that are below the filter,
  # and then extract the relevant data
  grouped_emf_peaks = purrr::map(seq(grouped_emf), function(ige){
    #message(ige)
    x = grouped_emf[[ige]]
    use_peaks = x$PeakID[1][[1]]
    peak_info = purrr::map(x$complete_EMF, ~ dplyr::filter(sample_assignments, PeakID %in% use_peaks, complete_EMF %in% .x))
    names(peak_info) = x$complete_EMF

    # the true e-value is the one that corresponds to the actual clique size, or what should be the maximum of the clique size
    grouped_evalues = purrr::map_dbl(x$complete_EMF, ~ dplyr::filter(evalue_clique_size, PeakID %in% use_peaks, complete_EMF %in% .x) %>% dplyr::slice(which.max(clique_size)) %>% dplyr::pull(e_value))
    grouped_clique_size = purrr::map_int(x$complete_EMF, ~ dplyr::filter(evalue_clique_size, PeakID %in% use_peaks, complete_EMF %in% .x) %>% dplyr::slice(which.max(clique_size)) %>% dplyr::pull(clique_size))

    keep_emf = grouped_evalues <= evalue_cutoff
    x = x[keep_emf, ]
    if (sum(keep_emf) > 0) {
      return(list(complete_EMF = x$complete_EMF,
                  peak_info = peak_info[keep_emf],
                  e_values = dplyr::mutate(x, e_value = grouped_evalues[keep_emf], Sample = sample_id, complete_EMF = complete_EMF, type = "primary"),
                  min_e_value = min(grouped_evalues[keep_emf]),
                  clique_size = grouped_clique_size[1],
                  PeakID = use_peaks,
                  Sample_Peak = paste0(peak_info[[1]]$Sample[1], "_", use_peaks)
      ))
    } else {
      return(NULL)
    }

  })

  # our initial list might have null's because they don't have any EMF assignments
  # with e-values below the cutoff
  gep_null = purrr::map_lgl(grouped_emf_peaks, is.null)
  grouped_emf_peaks = grouped_emf_peaks[!gep_null]
  names(grouped_emf_peaks) = paste0("GEMF_", seq_along(grouped_emf_peaks), ".", sample_id)

  # now we can use the isotopologue_EMF to find those things that have multiple adducts
  # to the same EMF, and then check if there is actually multiple evidences
  emf_by_emf = purrr::map_df(grouped_emf_peaks, function(x){
    purrr::map_dfr(x$peak_info, ~ data.frame(complete_EMF = .x$complete_EMF[1],
                                             adduct_IMF = .x$adduct_IMF[1],
                                             isotopologue_EMF = unique(dplyr::filter(.x, Type %in% "isotopologue_EMF") %>% dplyr::pull(Assignment_Data)), stringsAsFactors = FALSE))
  }) %>% split(., .$isotopologue_EMF)

  # can only have multiple evidences if there are multiple rows from above,
  # and H vs NH4 don't count as multiple evidence. This actually goes through
  # to find those instances
  n_adduct = purrr::map_dbl(emf_by_emf, nrow)
  emf_by_emf = emf_by_emf[n_adduct > 1]
  multi_evidence_emf = purrr::map_df(emf_by_emf, function(x){
    out_data = NULL
    if (length(unique(x$adduct_IMF)) > 1) {
      #emf_adducts = stringr::str_split_fixed(x$emf_Adduct, "\\.", 2)

      if (sum(x$adduct_IMF %in% c("1H1", "14N1,1H4")) != length(x$adduct_IMF)) {
        out_data = x
        out_data
      }
    }
    out_data
  })

  if (use_corroborating) {
    # extract the e-values for each one
    all_evalues = purrr::map_df(grouped_emf_peaks, ~ .x$e_values)

    gep_2_complete_emf = purrr::map2_dfr(grouped_emf_peaks, names(grouped_emf_peaks),
                                         function(.x, .y){
                                           #message(.y)
                                           data.frame(gep = .y, complete_EMF = .x$e_values$complete_EMF,
                                                      stringsAsFactors = FALSE)
                                         })
    gep_2_complete_emf = dplyr::filter(gep_2_complete_emf, complete_EMF %in% multi_evidence_emf$complete_EMF)
    unique_gep = unique(gep_2_complete_emf$gep)
    has_multi_evidence = which(names(grouped_emf_peaks) %in% unique_gep)

    # add in the e-values, with the note that they are from "corroborating" evidence,
    # so the actual assignments can be filtered out at the voting step.
    grouped_emf_peaks = purrr::map_at(grouped_emf_peaks, has_multi_evidence, function(gep){
      curr_evalues = gep$e_values

      multi_evalues = purrr::map_df(seq_len(nrow(curr_evalues)), function(curr_row){
        base_data = curr_evalues[curr_row, ]

        use_iso_emf = dplyr::filter(multi_evidence_emf, complete_EMF %in% base_data$complete_EMF) %>%
          dplyr::pull(isotopologue_EMF)
        emf_matches = dplyr::filter(multi_evidence_emf, !(complete_EMF %in% base_data$complete_EMF), isotopologue_EMF %in% use_iso_emf) %>%
          dplyr::pull(complete_EMF)

        if (length(emf_matches) > 0) {
          other_data = dplyr::filter(all_evalues, complete_EMF %in% emf_matches) %>%
            dplyr::mutate(complete_EMF = base_data$complete_EMF, type = "corroborating")
          base_data = rbind(base_data, other_data)
        }
        base_data
      })
      gep$e_values = multi_evalues
      gep

    })
  }


  return(list(grouped_emf = grouped_emf_peaks, multi_evidence = multi_evidence_emf))
  # we should add the peaks to GEMF mapping here, so we can do what is noted on the last line.
}

#' create sudo elemental molecular formulas
#'
#' Given the grouped EMFs in a sample to EMFs, come up with `sudo` EMFs, that is groups
#' of shared EMFs within and across samples, based on shared EMFs.
#'
#' @param gemf_2_emf list of grouped EMFs to EMFs
#'
#' @details This is achieved by repeated set intersection / set union cycles, where
#'  they are repeated until the set of things is consistent between cycles, and there
#'  are no entries to pull from.
#'
#' @return list
#' @export
create_sudo_emfs = function(gemf_2_emf){
  sudo_emf_list = vector("list", length(unique(gemf_2_emf$grouped_emf)))

  i_sudo = 1
  while (nrow(gemf_2_emf) > 0) {
    #message(i_sudo)
    tmp_emf_mapping = dplyr::filter(gemf_2_emf, complete_EMF %in% gemf_2_emf$complete_EMF[1])

    emf_iter = dplyr::left_join(tmp_emf_mapping, gemf_2_emf, by = "grouped_emf") %>%
      dplyr::transmute(grouped_emf = grouped_emf, complete_EMF = complete_EMF.y) %>% unique()

    adduct_iter = dplyr::left_join(emf_iter, gemf_2_emf, by = "complete_EMF") %>%
      dplyr::transmute(grouped_emf = grouped_emf.y, complete_EMF = complete_EMF) %>% unique()

    while (nrow(adduct_iter) != nrow(emf_iter)) {
      emf_iter = dplyr::left_join(adduct_iter, gemf_2_emf, by = "grouped_emf") %>%
        dplyr::transmute(grouped_emf = grouped_emf, complete_EMF = complete_EMF.y) %>% unique()

      adduct_iter = dplyr::left_join(emf_iter, gemf_2_emf, by = "complete_EMF") %>%
        dplyr::transmute(grouped_emf = grouped_emf.y, complete_EMF = complete_EMF) %>% unique()
    }

    sudo_emf_list[[i_sudo]] = adduct_iter
    gemf_2_emf = dplyr::filter(gemf_2_emf, !(grouped_emf %in% adduct_iter$grouped_emf))
    i_sudo = i_sudo + 1
  }

  keep_sudo = purrr::map_lgl(sudo_emf_list, ~ !is.null(.x))
  sudo_emf_list = sudo_emf_list[keep_sudo]
  names(sudo_emf_list) = paste0("SEMF.", seq(1, length(sudo_emf_list)))
  sudo_emf_list

}


#' match peaks to IMFs by M/Z
#'
#' Given a data.frame of IMFs to peaks for a given EMF, attempts to match unknown peaks
#' to IMFs based on differences in average M/Z to the unknown peaks.
#'
#' @param imf_2_peak data.frame with IMF information (see details)
#' @param unknown_peaks data.frame of peak information
#' @param peak_mz data.frame with `Sample_Peak` and `Value`, where `Value` is `ObservedMZ`
#'
#' @return data.frame
#' @export
match_imf_by_mz = function(imf_2_peak, unknown_peaks, peak_mz){
  split_peaks_imf = split(imf_2_peak$Sample_Peak, imf_2_peak$complete_IMF)
  mean_mz = purrr::map2_dfr(split_peaks_imf, names(split_peaks_imf),
                            function(.x, .y){
                              use_mz = peak_mz[peak_mz$Sample_Peak %in% .x, "Value"]
                              data.frame(complete_IMF = .y,
                                         mean = mean(use_mz),
                                         sd = 2*sd(use_mz),
                                         stringsAsFactors = FALSE)
                            })
  mean_mz = dplyr::mutate(mean_mz, ppm = mean / 1e6)
  mean_mz = dplyr::mutate(mean_mz, sd2 = dplyr::case_when(sd > ppm ~ ppm,
                                                          sd <= ppm ~ sd))

  peak_nap = unique(imf_2_peak[, c("complete_IMF", "NAP")])
  mean_mz = dplyr::left_join(mean_mz, peak_nap, by = "complete_IMF") %>%
    dplyr::arrange(dplyr::desc(NAP)) %>% dplyr::mutate(seq = seq(1, dplyr::n()))

  just_peaks = as.character(unknown_peaks$Peaks)
  # this also needs to have the peaks sorted by their NAP, and match in order
  # of NAP. When the next one doesn't match, stop!
  match_mz = purrr::map_df(just_peaks, function(test_peak){
    test_mz = dplyr::filter(peak_mz, Sample_Peak %in% test_peak)
    tmp_mean = mean_mz %>% dplyr::mutate(diff = abs(mean - test_mz$Value)) %>%
      dplyr::filter(diff <= sd2)

    if (nrow(tmp_mean) > 0) {
      tmp_mean = dplyr::slice(tmp_mean, which.min(diff))
      match_imf = data.frame(complete_IMF = tmp_mean$complete_IMF,
                             PeakID = test_mz$PeakID,
                             Sample = test_mz$Sample,
                             Sample_Peak = test_peak,
                             complete_EMF = imf_2_peak$complete_EMF[1],
                             seq = tmp_mean$seq[1],
                             stringsAsFactors = FALSE)
    } else {
      match_imf = data.frame(complete_IMF = as.character(NA),
                             PeakID = NA,
                             Sample = as.character(NA),
                             Sample_peak = as.character(NA),
                             complete_EMF = as.character(NA),
                             seq = as.integer(NA),
                             stringsAsFactors = FALSE)
    }
    match_imf
  })

  match_mz = match_mz[!is.na(match_mz$complete_IMF), ]

  null_evalue = data.frame(complete_EMF = imf_2_peak$complete_EMF[1], e_value = NA, stringsAsFactors = FALSE)
  if (nrow(match_mz) > 0) {
    match_mz = dplyr::arrange(match_mz, seq)
    seq_match = all(match_mz$seq == seq_len(nrow(match_mz)))
    unique_match = length(unique(match_mz$complete_IMF)) == nrow(match_mz)

    if (!(seq_match && unique_match)) {
      #print(match_mz)
      match_mz = match_mz[1,]
      match_mz[1, ] = NA
    }
  } else {
    match_mz = match_mz[1,]
    match_mz[1, ] = NA
  }

  out_evalues = null_evalue
  out_evalues$info = list(match_mz)
  out_evalues$Sample_Peak = list(match_mz$Sample_Peak)
  out_evalues$Sample = match_mz$Sample[1]
  out_evalues$grouped_emf = as.character(unknown_peaks$grouped_emf[1])

  out_evalues
}

#' match peaks to IMFs by difference
#'
#' Given a data.frame of IMFs to peaks for a given EMF, attempts to match unknown peaks
#' to IMFs based on differences in the provided values to the unknown peaks.
#'
#' @param imf_2_peak data.frame with IMF information (see details)
#' @param unknown_peaks data.frame of peak information
#' @param scan_locations list of locations for each peak across scans
#' @param peak_location data.frame with `Sample_Peak` and `Value`, where `Value` is `ObservedFrequency`
#' @param difference_cutoff the maximum difference an unknown peak can have to any of the known peaks
#'
#' @return data.frame
#' @export
match_imf_by_difference = function(imf_2_peak, unknown_peaks, scan_locations, peak_location, difference_cutoff){
  split_peaks_imf = split(imf_2_peak$Sample_Peak, imf_2_peak$complete_IMF)
  split_peaks_location = purrr::map(split_peaks_imf, function(in_imf){
    peak_location[peak_location$Sample_Peak %in% in_imf, "Value"]
  })

  # if we didn't supply a difference cutoff, then use 1ppm. Should work in
  # either space.
  if (is.na(difference_cutoff)) {
    "!DEBUG using 1 ppm or 1 frequency cutoff"
    if ("ObservedMZ" %in% names(difference_cutoff)) {
      max_loc = max(unlist(split_peaks_location))
      difference_cutoff = max_loc / 1e6
    } else if ("ObservedFrequency" %in% names(difference_cutoff)) {
      difference_cutoff = 1
    }

  }

  mean_location = purrr::imap_dfr(split_peaks_location, function(peak_locs, in_imf){
    data.frame(location = mean(peak_locs),
               peak_sd = 2*sd(peak_locs),
               n_peak = length(peak_locs),
               imf = in_imf,
               stringsAsFactors = FALSE)
  })

  mean_location = dplyr::mutate(mean_location,
      peak_sd2 = dplyr::case_when(is.na(peak_sd) ~ difference_cutoff,
                                  peak_sd > difference_cutoff ~ difference_cutoff,
                                  peak_sd <= difference_cutoff ~ peak_sd))

  if (!is.null(scan_locations)) {
    "!DEBUG using scan level data"
    imf_scan_location_sd = purrr::imap_dfr(split_peaks_imf, function(imf_peaks, in_imf){
      tmp_scan_loc = unlist(scan_locations[imf_peaks])
      data.frame(scan_sd = sd(tmp_scan_loc, na.rm = TRUE) * 2,
                 imf = in_imf,
                 stringsAsFactors = FALSE)
    })

    mean_location = dplyr::left_join(mean_location, imf_scan_location_sd, by = "imf")

    if (max(mean_location$n_peak) > 2) {
      if (min(mean_location$n_peak) < 2) {
        other_ranges = dplyr::filter(mean_location, n_peak > 2) %>% dplyr::pull(scan_sd) %>%
          max()
        low_loc = which(mean_location$n_peak < 2)
        for (iloc in low_loc) {
          if (mean_location[iloc, "scan_sd"] < (0.9 * other_ranges)) {
            mean_location[iloc, "scan_sd"] = difference_cutoff
          }
        }

      }
    }

    mean_location = dplyr::mutate(mean_location,
      use_sd = purrr::pmap_dbl(mean_location[, c("scan_sd", "peak_sd2")],
                               ~min(c(...), na.rm = TRUE)))

  } else {
    "!DEBUG no scan level data used"
    mean_location$use_sd = mean_location$peak_sd2
  }

  peak_nap = unique(imf_2_peak[, c("complete_IMF", "NAP")]) %>%
    dplyr::arrange(dplyr::desc(NAP)) %>% dplyr::mutate(seq = seq(1, dplyr::n()))
  mean_location = dplyr::left_join(mean_location, peak_nap, by = c("imf" = "complete_IMF"))
  split_peaks_location = split_peaks_location[peak_nap$complete_IMF]

  just_peaks = as.character(unknown_peaks$Peaks)
  # this also needs to have the peaks sorted by their NAP, and match in order
  # of NAP. When the next one doesn't match, stop!
  match_location = purrr::map_df(just_peaks, function(test_peak){
    test_location = dplyr::filter(peak_location, Sample_Peak %in% test_peak)
    mean_diffs = mean_location
    mean_diffs$diff = abs(mean_location$location - test_location$Value)

    mean_diffs = dplyr::filter(mean_diffs, diff <= use_sd)
    if (nrow(mean_diffs) > 0) {
      tmp_location = dplyr::slice(mean_diffs, which.min(diff))
      match_imf = data.frame(complete_IMF = tmp_location$imf,
                             PeakID = test_location$PeakID,
                             Sample = test_location$Sample,
                             Sample_Peak = test_peak,
                             complete_EMF = imf_2_peak$complete_EMF[1],
                             NAP = tmp_location$NAP,
                             seq = tmp_location$seq,
                             stringsAsFactors = FALSE)
    } else {
      match_imf = data.frame(complete_IMF = as.character(NA),
                             PeakID = NA,
                             Sample = as.character(NA),
                             Sample_Peak = as.character(NA),
                             complete_EMF = as.character(NA),
                             NAP = as.double(NA),
                             seq = as.integer(NA),
                             stringsAsFactors = FALSE)
    }
    match_imf
  })

  match_location = match_location[!is.na(match_location$complete_IMF), ]

  null_evalue = data.frame(complete_EMF = imf_2_peak$complete_EMF[1], e_value = NA, stringsAsFactors = FALSE)
  if (nrow(match_location) > 0) {
    match_location = dplyr::arrange(match_location, seq)
    seq_match = all(match_location$seq == seq_len(nrow(match_location)))
    unique_match = length(unique(match_location$complete_IMF)) == nrow(match_location)

    if (!(seq_match && unique_match)) {
      #print(match_location)
      match_location = match_location[1,]
      match_location[1, ] = NA
    }
  } else {
    match_location = match_location[1,]
    match_location[1, ] = NA
  }

  out_evalues = null_evalue
  out_evalues$info = list(match_location)
  out_evalues$Sample_Peak = list(match_location$Sample_Peak)
  out_evalues$Sample = match_location$Sample[1]
  out_evalues$grouped_emf = as.character(unknown_peaks$grouped_emf[1])

  out_evalues
}


#' choose most likely elemental molecular formula's
#'
#' Given a list of grouped EMFs (assumed to be from a single sudo EMF), votes on the most likely
#' EMF based on the sum of `1 - e_value` across the grouped EMFs. If a grouped EMF does not
#' have a matching EMF, then attempts to match peaks to the voted EMFs based on M/Z and matching
#' to EMF peaks in natural abundance order.
#'
#' @param grouped_emfs list of grouped EMFs
#' @param peak_location data.frame of peak locations
#' @param difference_cutoff how close do peaks need to be to each other
#' @param keep_ratio how close to the maximum voted EMF to keep other things? (default is 0.9)
#'
#' @return data.frame
#' @export
choose_emf = function(grouped_emfs, scan_level_location, peak_location, difference_cutoff, keep_ratio = 0.9){
  grouped_evalues = purrr::map_df(grouped_emfs, ~ .x$e_values) %>% dplyr::mutate(information = 1 - e_value)

  use_peaks = unique(unlist(purrr::map(grouped_emfs, "Sample_Peak")))
  use_location = dplyr::filter(peak_location, Sample_Peak %in% use_peaks)

  if (!is.null(scan_level_location)) {
    use_scans = scan_level_location[use_peaks]
  } else {
    use_scans = NULL
  }

  if (!any(is.na(difference_cutoff)) && is.data.frame(difference_cutoff)) {
    max_location = max(use_location$Value)

    difference_loc = which.min(abs(max_location - difference_cutoff$Index))
    use_difference_cutoff = difference_cutoff[difference_loc, "Value"]
  } else {
    use_difference_cutoff = difference_cutoff
  }


  # sum 1-evalue across the samples for each EMF, and keep those things that are
  # within X% of the highest sum
  emf_votes = dplyr::group_by(grouped_evalues, complete_EMF) %>% dplyr::summarise(sum_information = sum(information)) %>%
    dplyr::mutate(max_ratio = sum_information / max(sum_information))
  keep_emf = dplyr::filter(emf_votes, max_ratio >= keep_ratio)

  null_evalue = data.frame(complete_EMF = NA, e_value = NA, stringsAsFactors = FALSE)
  gemf_emf = purrr::map2_dfr(grouped_emfs, names(grouped_emfs), function(.x, .y){
    sample_from_grouped_emf = stringr::str_split_fixed(.y, "\\.", 2)[2]
    e_values = dplyr::filter(.x$e_values, complete_EMF %in% keep_emf$complete_EMF, type %in% "primary") %>% dplyr::select(complete_EMF, e_value)

    if (nrow(e_values) > 0) {
      e_values$info = vector(mode = "list", length = nrow(e_values))
      e_values$Sample_Peak = vector(mode = "list", length = nrow(e_values))
      use_info = .x$peak_info[e_values$complete_EMF]
      e_values$info = purrr::map(use_info, function(in_info){
        in_info %>% dplyr::filter(Type %in% "NAP") %>% dplyr::mutate(NAP = as.numeric(Assignment_Data)) %>%
          dplyr::select(complete_IMF, PeakID, Sample, Sample_Peak, complete_EMF, NAP) %>%
          unique()
      })

      e_values$Sample_Peak = list(.x$Sample_Peak)
      e_values$Sample = as.character(sample_from_grouped_emf)
      e_values$grouped_emf = as.character(.y)
    } else {
      e_values = null_evalue
      e_values$info = vector(mode = "list", length = nrow(e_values))
      e_values$Sample_Peak = list(.x$Sample_Peak)
      e_values$Sample = as.character(sample_from_grouped_emf)
      e_values$grouped_emf = as.character(.y)
    }
    e_values
  })
  na_evalues = purrr::map_lgl(gemf_emf$e_value, is.na)

  # if (sum(na_evalues) > 0) {
  #   message("Found one!")
  # }

  have_imf_peaks = unique(unlist(purrr::map(gemf_emf$Sample_Peak[!na_evalues], ~ .x), use.names = FALSE))

  if (sum(na_evalues) > 0) {
    trimmed_gemf_emf = gemf_emf[!na_evalues, ]
    missing_imf = purrr::map(which(na_evalues), function(.x){
      tmp_df = data.frame(Peaks = as.character(gemf_emf$Sample_Peak[[.x]]), Sample = gemf_emf$Sample[.x],
                          grouped_emf = gemf_emf$grouped_emf[.x], stringsAsFactors = FALSE)
      if (all(tmp_df$Peaks %in% have_imf_peaks)) {
        return(NULL)
      } else {
        return(tmp_df)
      }
    })
    missing_imf = missing_imf[!purrr::map_lgl(missing_imf, is.null)]
    out_gemf_emf = trimmed_gemf_emf
    if (length(missing_imf) > 0) {

      has_imf = unique(purrr::map_df(trimmed_gemf_emf$info, ~ .x))

      imf_by_emf = split(has_imf, has_imf$complete_EMF)

      imf_matches = purrr::map_df(purrr::cross2(imf_by_emf, missing_imf), function(x){
        match_imf_by_difference(x[[1]], x[[2]], use_scans, use_location, use_difference_cutoff)
      })
      imf_matches = imf_matches[!is.na(imf_matches$Sample), ]

      if (nrow(imf_matches) > 0) {
        out_gemf_emf = rbind(out_gemf_emf, imf_matches)
      }
    }
  } else {
    out_gemf_emf = gemf_emf
  }

  # ideally, right here we should go through and confirm that the the combination
  # of EMF and Sample_Peak is unique, and remove any duplicate / subsumed entries.
  # this seems like a hard problem to solve right now, so leaving it for later.

  # finally, go through each EMF and the associated peaks, and confirm that *all*
  # of the peaks for each of the IMFs are within the difference_cutoff, and
  # remove any that aren't
  if (!is.na(use_difference_cutoff)) {
    "!DEBUG checking imfs"
    checked_emfs = purrr::map(split(out_gemf_emf, out_gemf_emf$complete_EMF), check_emf, use_location, use_difference_cutoff)
    returned_emfs = do.call(rbind, checked_emfs)
  } else {
    "!DEBUG No checking imfs"
    returned_emfs = out_gemf_emf
  }

  returned_emfs
}

check_emf = function(in_emf, peak_location, difference_cutoff){
  emf_info = purrr::map_df(in_emf$info, ~ .x)

  peaks_by_imf = split(emf_info$Sample_Peak, emf_info$complete_IMF)

  keep_peaks_by_imf = purrr::map(peaks_by_imf, function(in_peaks){
    in_frequencies = dplyr::filter(peak_location, Sample_Peak %in% in_peaks)

    if (nrow(in_frequencies) == 1) {
      return(in_peaks)
    } else {
      median_diffs = purrr::map_dbl(seq(1, nrow(in_frequencies)), function(in_row){
        median(abs(in_frequencies[in_row, "Value"] - in_frequencies[-in_row, "Value"]))
      })
      out_peaks = in_frequencies[median_diffs <= difference_cutoff, "Sample_Peak"]
      return(out_peaks)
    }
  })

  all_keep = unique(unlist(keep_peaks_by_imf))

  info_by_sample = split(emf_info, emf_info$Sample) %>%
    purrr::map_df(., function(in_sample){
      in_sample = in_sample[order(in_sample$NAP, decreasing = TRUE), ]
      in_sample$order = seq(1, nrow(in_sample))

      in_sample = dplyr::filter(in_sample, Sample_Peak %in% all_keep)

      if (nrow(in_sample) > 0) {
        tmp_keep = rep(FALSE, nrow(in_sample))
        for (i_peak in seq_len(nrow(in_sample))) {
          if (in_sample[i_peak, "order"] == 1) {
            tmp_keep[i_peak] = TRUE
          } else {
            if ((in_sample[i_peak, "order"] - 1) %in% in_sample[tmp_keep, "order"]) {
              tmp_keep[i_peak] = TRUE
            }
          }
        }
        in_sample = in_sample[tmp_keep, ]
      }
      in_sample
    })

  in_emf$info = purrr::map(in_emf$info, function(in_info){
    dplyr::filter(in_info, Sample_Peak %in% info_by_sample$Sample_Peak)
  })
  in_emf$Sample_Peak = purrr::map(in_emf$Sample_Peak, function(in_peaks){
    base::intersect(in_peaks, info_by_sample$Sample_Peak)
  })

  nonzero_rows = purrr::map_lgl(in_emf$info, ~ nrow(.x) > 0)
  in_emf = in_emf[nonzero_rows, ]
  in_emf
}

#' remove duplicate peaks across sudo EMFs
#'
#' Given a set of sudo EMFs, detects duplicate peaks and removes them, and carries
#' out EMF voting on what is left. Should only be used after an initial round of
#' merging.
#'
#' @param chosen_emfs merged EMFs
#' @param all_gemfs the grouped_EMFs
#' @param peak_location the data.frame of peak_location
#' @param difference_cutoff tolerance for matched peaks for `choose_emf`
#' @param keep_ratio for `choose_emf`
#'
#' @return list of sudo EMFs after voting
#' @export
remove_duplicates_across_semfs = function(chosen_emfs, all_gemfs, scan_level_location, peak_location, difference_cutoff, keep_ratio = 0.9){
  peak_2_voted_emf = purrr::map2_df(chosen_emfs, names(chosen_emfs), function(.x, .y){
    #message(.y)
    data.frame(Sample_Peak = unique(unlist(.x$Sample_Peak, use.names = FALSE)),
               semf = .y,
               stringsAsFactors = FALSE)
  })

  dup_peaks = unique(peak_2_voted_emf$Sample_Peak[duplicated(peak_2_voted_emf$Sample_Peak)])
  has_dup_semf = dplyr::filter(peak_2_voted_emf, Sample_Peak %in% dup_peaks)

  has_dup_emf = which(names(chosen_emfs) %in% unique(has_dup_semf$semf))
  non_dup_emf = which(!(names(chosen_emfs) %in% unique(has_dup_semf$semf)))

  non_dup_chosen = chosen_emfs[non_dup_emf]
  potential_dup_chosen = chosen_emfs[has_dup_emf]

  chosen_emfs_dedup = internal_map$map_function(potential_dup_chosen, function(in_semf){
    peak_2_gemf = purrr::map_df(seq(1, nrow(in_semf)), function(in_row){
      data.frame(Sample_Peak = in_semf[in_row, "Sample_Peak"][[1]],
                 grouped_emf = in_semf[in_row, "grouped_emf"],
                 stringsAsFactors = FALSE)
    })
    bad_gemfs = dplyr::filter(peak_2_gemf, Sample_Peak %in% dup_peaks)
    good_gemfs = dplyr::filter(peak_2_gemf, !(grouped_emf %in% bad_gemfs$grouped_emf))

    if (nrow(good_gemfs) > 0) {
      return(choose_emf(all_gemfs[unique(good_gemfs$grouped_emf)], scan_level_location, peak_location, difference_cutoff, keep_ratio))
    } else {
      return(NULL)
    }

  })

  dedup_null = purrr::map_lgl(chosen_emfs_dedup, is.null)
  out_emfs = c(non_dup_chosen, chosen_emfs_dedup[!dedup_null])
  names(out_emfs) = paste0("SEMF.", seq(1, length(out_emfs)))
  out_emfs
}


#' merge sudo EMFs with duplicate peaks
#'
#' After choosing EMFs, there may be sudo EMFs that share peaks. This is not kosher, there
#' has to be a 1-1 mapping of peaks to sudo EMFs. This function merges those that share
#' > 50% of peaks from either sudo EMF, and carries out voting again.
#' After merging, it checks for duplicates again, and removes those peaks that still
#' span multiple sudo EMFs, and then votes on them again.
#'
#' @param chosen_emfs sudo EMFs from `choose_emf`
#' @param all_gemfs the grouped_EMFs
#' @param peak_location frequency information for each and every peak
#' @param difference_cutoff how close do peaks need to be to each other
#' @param keep_ratio the ratio used to keep other highly voted EMFs
#'
#' @return list
#' @export
merge_duplicate_semfs = function(chosen_emfs, all_gemfs, scan_level_location, peak_location, difference_cutoff, keep_ratio = 0.9){
  peak_2_voted_emf = purrr::map2_df(chosen_emfs, names(chosen_emfs), function(.x, .y){
    all_peaks = unique(unlist(purrr::map(.x$Sample_Peak, ~ .x), use.names = FALSE))
    data.frame(Sample_Peak = all_peaks, semf = .y, stringsAsFactors = FALSE)
  })


  dup_peaks = dplyr::filter(peak_2_voted_emf, Sample_Peak %in% (unique(peak_2_voted_emf$Sample_Peak[duplicated(peak_2_voted_emf$Sample_Peak)]))) %>%
    dplyr::arrange(Sample_Peak)

  use_semfs = unique(dup_peaks$semf)

  index_semfs = dup_semfs = purrr::map(use_semfs, ~ .x)
  names(index_semfs) = names(dup_semfs) = use_semfs
  index_semfs_peaks = purrr::map(chosen_emfs[use_semfs], ~ unique(unlist(.x$Sample_Peak, use.names = FALSE)))
  index_semfs_gemfs = purrr::map(chosen_emfs[use_semfs], ~ unique(unlist(.x$grouped_emf, use.names = FALSE)))

  dup_semfs_peaks = index_semfs_peaks
  dup_semfs_gemfs = index_semfs_gemfs

  for (isemf in use_semfs) {
    use_i_semfs = index_semfs[[isemf]]
    i_peaks = unique(unlist(index_semfs_peaks[use_i_semfs], use.names = FALSE))
    n_i_peaks = length(i_peaks)
    for (jsemf in use_semfs[!(use_semfs %in% isemf)]) {
      use_j_semfs = index_semfs[[jsemf]]
      j_peaks = unique(unlist(index_semfs_peaks[use_j_semfs], use.names = FALSE))
      n_j_peaks = length(j_peaks)
      intersect_ratio = length(intersect(i_peaks, j_peaks)) / c(n_i_peaks, n_j_peaks)

      if (max(intersect_ratio) >= 0.5) {
        all_semfs = unique(c(use_i_semfs, use_j_semfs))
        index_semfs[[isemf]] = index_semfs[[jsemf]] = all_semfs
        dup_semfs_peaks[all_semfs] = list(unique(unlist(index_semfs_peaks[all_semfs], use.names = FALSE)))
        dup_semfs_gemfs[all_semfs] = list(unique(unlist(index_semfs_gemfs[all_semfs], use.names = FALSE)))
      }
    }
  }
  keep_semfs = !duplicated(index_semfs)
  merged_gemfs = purrr::map(dup_semfs_gemfs[keep_semfs], ~ data.frame(grouped_emf = .x, stringsAsFactors = FALSE))
  dup_semfs = use_semfs

  chosen_nondup = chosen_emfs[!(names(chosen_emfs) %in% use_semfs)]
  chosen_duplicates = internal_map$map_function(merged_gemfs, function(.x){
    choose_emf(all_gemfs[unique(.x$grouped_emf)], scan_level_location, peak_location, difference_cutoff, keep_ratio)
  })

  all_chosen_emfs = c(chosen_nondup, chosen_duplicates)
  names(all_chosen_emfs) = paste0("SEMF.", seq(1, length(all_chosen_emfs)))

  remove_duplicates_across_semfs(all_chosen_emfs, all_gemfs, scan_level_location, peak_location, difference_cutoff, keep_ratio)

}

extract_emfs = function(chosen_emfs){
  all_samples = unique(unlist(purrr::map(chosen_emfs, ~ .x$Sample)))

  null_sample_matrix = matrix(as.character(NA), nrow = 1, ncol = length(all_samples))
  colnames(null_sample_matrix) = all_samples

  #all_emfs = internal_map$map_function(chosen_emfs, function(use_emf){
  all_emfs = purrr::map(seq(1, length(chosen_emfs)), function(i_emf){
    #message(i_emf)
    use_emf = chosen_emfs[[i_emf]]
    peak_info = purrr::map_df(use_emf$info, ~ .x)

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

  semf_emf = purrr::map2_df(all_emfs, names(all_emfs), function(.x, .y){
    data.frame(semf = .y, emf = names(.x), stringsAsFactors = FALSE)
  })

  if (sum(duplicated(semf_emf$emf)) > 0) {
    warning("Duplicate sudo EMF to complete EMF mappings detected!")
  }

  all_emfs
}

#' extract IMF or EMF level data
#'
#' Extract either IMF or EMF level data from the extracted EMFs. Which is controlled by the `by` argument.
#'
#' @param emfs list of EMFs
#' @param by `EMF` (default) or `IMF`
#'
#' @details For `EMF`, will be a list corresponding to each EMF from the sudo EMFs, for `IMF`, will be
#'  matrices and a data.frame.
#'
#' @return list with `height`, `mz`, and `info`
#' @export
#'
extract_imf_emf_data = function(emfs, by = "EMF"){
  all_compared = internal_map$map_function(emfs, compare_extract_emfs, by = by)

  if (is.null(names(all_compared))) {
    names(all_compared) = paste0("SEMF.", seq_along(all_compared))
  }

  if ("EMF" %in% by) {
    all_height = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_height = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, "_", emf_id, "_", names(in_emf$info))
        out_height = in_emf$height
        rownames(out_height) = peak_id
        out_height
      })
      names(semf_height) = paste0(semf_id, "_", names(in_semf))
      semf_height

    }) %>% purrr::flatten()

    all_mz = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_mz = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, "_", emf_id, "_", names(in_emf$info))
        out_mz = in_emf$mz
        rownames(out_mz) = peak_id
        out_mz
      })
      names(semf_mz) = paste0(semf_id, "_", names(in_semf))
      semf_mz

    }) %>% purrr::flatten()

    all_info = purrr::map2(all_compared, names(all_compared), function(in_semf, semf_id){
      semf_info = purrr::map2(in_semf, names(in_semf), function(in_emf, emf_id){
        peak_id = paste0(semf_id, "_", emf_id, "_", names(in_emf$info))
        info_df = purrr::map2_dfr(in_emf$info, peak_id, function(in_info, pid){
          in_info$PeakID = pid
          in_info$sudo_EMF = semf_id
          in_info$sudo_group_EMF = paste0(emf_id)
          in_info
        })
      })
      names(semf_info) = paste0(semf_id, "_", names(in_semf))
      semf_info
    }) %>% purrr::flatten()

  } else {
    all_height = purrr::map2(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,"_", names(.x$info))
      out_height = .x$height
      n_height = nrow(out_height)
      message(paste0(.y, " ", n_height))
      rownames(out_height) = peak_id
      out_height
    }) %>% do.call(rbind, .)
    all_mz = purrr::map2(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,"_", names(.x$info))
      out_mz = .x$mz
      rownames(out_mz) = peak_id
      out_mz
    }) %>% do.call(rbind, .)
    all_info = purrr::map2_dfr(all_compared, names(all_compared), function(.x, .y){
      peak_id = paste0(.y,"_", names(.x$info))
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
  return(list(height = all_height,
              mz = all_mz,
              info = all_info))
}

compare_extract_emfs = function(emf, by = "EMF"){

  get_imfs = function(emf_group){
    split_groups = split(emf_group, emf_group$group_imfs)
    names(split_groups) = paste0("IMF.", seq_along(split_groups))
    single_index = purrr::map_df(split_groups, ~ .x[1, ])

    out_height = compare_heights[single_index$row_index, ]
    out_mz = compare_mz[single_index$row_index, ]

    out_info = purrr::map(split_groups, ~ compare_info[.x$row_index, ])

    list(height = out_height, mz = out_mz, info = out_info)
  }

  compare_heights = purrr::map(emf, ~ .x$height) %>% do.call(rbind, .)
  compare_info = purrr::map_df(emf, ~ .x$peak_info)
  compare_mz = purrr::map(emf, ~ .x$mz) %>% do.call(rbind, .)

  group_imfs = vector("integer", length = nrow(compare_heights))

  i_group = 1
  while (sum(group_imfs == 0) > 0) {
    zero_locs = which(group_imfs == 0)
    master_loc = which(group_imfs == 0)[1]

    for (iloc in zero_locs) {
      if (isTRUE(all.equal(compare_heights[master_loc, ], compare_heights[iloc, ]))) {
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
    names(split_emfs) = paste0("EMF.", seq_along(split_emfs))

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
#' @param emf_heights list of EMF heights
#'
#' @export
#' @return matrix
#'
extract_emf_propensity = function(emf_heights){
  presence_matrix = purrr::map(emf_heights, ~ as.integer(colSums(.x) > 0)) %>% do.call(rbind, .)
  rownames(presence_matrix) = names(emf_heights)
  colnames(presence_matrix) = colnames(emf_heights[[1]])
  presence_matrix
}
