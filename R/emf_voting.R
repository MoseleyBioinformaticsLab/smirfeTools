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
  sample_assignments = dplyr::mutate(sample_assignments, emf_Adduct = paste0(complete_EMF, ".", adduct_IMF))
  e_values = dplyr::filter(sample_assignments, Type %in% "e_value") %>% dplyr::mutate(e_value = as.numeric(Assignment_Data), imf_peak = paste0(complete_IMF, "_", adduct_IMF, "_", PeakID))
  clique_size = dplyr::filter(sample_assignments, Type %in% "clique_size") %>% dplyr::mutate(clique_size = as.integer(Assignment_Data), imf_peak = paste0(complete_IMF, "_", adduct_IMF, "_", PeakID))
  evalue_clique_size = dplyr::left_join(e_values, clique_size[, c("imf_peak", "clique_size")], by = "imf_peak")

  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$complete_EMF) %>%
    purrr::map2_dfr(., names(.), function(.x, .y){
      peaks = unique(.x)
      tmp_frame = data.frame(complete_EMF = .y, PeakID_chr = paste(peaks, collapse = ","), stringsAsFactors = FALSE)
      tmp_frame$PeakID = list(peaks)
      tmp_frame
    })

  grouped_emf = split(peaks_by_emf, peaks_by_emf$PeakID_chr)

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

  gep_null = purrr::map_lgl(grouped_emf_peaks, is.null)
  grouped_emf_peaks = grouped_emf_peaks[!gep_null]
  names(grouped_emf_peaks) = paste0("GEMF_", seq_along(grouped_emf_peaks), ".", sample_id)

  emf_by_emf = purrr::map_df(grouped_emf_peaks, function(x){
    purrr::map_dfr(x$peak_info, ~ data.frame(complete_EMF = .x$complete_EMF[1],
                                             adduct_IMF = .x$adduct_IMF[1],
                                             isotopologue_EMF = unique(dplyr::filter(.x, Type %in% "isotopologue_EMF") %>% dplyr::pull(Assignment_Data)), stringsAsFactors = FALSE))
  }) %>% split(., .$isotopologue_EMF)

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

#' choose most likely elemental molecular formula's
#'
#' Given a list of grouped EMFs (assumed to be from a single sudo EMF), votes on the most likely
#' EMF based on the sum of `1 - e_value` across the grouped EMFs. If a grouped EMF does not
#' have a matching EMF, then attempts to match peaks to the voted EMFs based on M/Z and matching
#' to EMF peaks in natural abundance order.
#'
#' @param grouped_emfs list of grouped EMFs
#' @param peak_mz data.frame of peak M/Zs
#' @param keep_ratio how close to the maximum voted EMF to keep other things? (default is 0.9)
#'
#' @return data.frame
#' @export
choose_emf = function(grouped_emfs, peak_mz, keep_ratio = 0.9){
  grouped_evalues = purrr::map_df(grouped_emfs, ~ .x$e_values) %>% dplyr::mutate(information = 1 - e_value)

  # I think the way to incorporate multi-adduct evidence is to add some small constant value to each
  # of the samples where that EMF also had corresponding evidence. Either that, or do an initial vote
  # with loose ratio cutoff (say 0.5), then filter for corresponding in each sample, and do the vote
  # again, with something added to boost the corresponding.

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

  have_imf_peaks = unique(unlist(purrr::map(gemf_emf$Sample_Peaks[!na_evalues], ~ .x), use.names = FALSE))

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
        match_imf_by_mz(x[[1]], x[[2]], peak_mz)
      })
      imf_matches = imf_matches[!is.na(imf_matches$sample), ]

      if (nrow(imf_matches) > 0) {
        out_gemf_emf = rbind(out_gemf_emf, imf_matches)
      }
    }
  } else {
    out_gemf_emf = gemf_emf
  }

  # if (length(setdiff(unique(gemf_emf$sample), unique(out_gemf_emf$sample))) > 0) {
  #   message("sample doesn't have IMFs")
  # }

  out_gemf_emf

}

#' remove duplicate peaks across sudo EMFs
#'
#' Given a set of sudo EMFs, detects duplicate peaks and removes them, and carries
#' out EMF voting on what is left. Should only be used after an initial round of
#' merging.
#'
#' @param chosen_emfs merged EMFs
#' @param all_gemfs the grouped_EMFs
#' @param peak_mz the data.frame of peak_mz
#' @param keep_ratio for `choose_emf`
#'
#' @return list of sudo EMFs after voting
#' @export
remove_duplicates_across_semfs = function(chosen_emfs, all_gemfs, peak_mz, keep_ratio = 0.9){
  peak_2_voted_emf = purrr::map2_df(chosen_emfs, names(chosen_emfs), function(.x, .y){
    #message(.y)
    data.frame(Sample_Peak = unique(unlist(.x$Sample_Peak, use.names = FALSE)),
               semf = .y,
               stringsAsFactors = FALSE)
  })

  dup_peaks = unique(peak_2_voted_emf$Sample_Peak[duplicated(peak_2_voted_emf$Sample_Peak)])
  has_dup_semf = dplyr::filter(peak_2_voted_emf, Sample_Peak %in% dup_peaks)

  has_dup_emf = which(names(chosen_emfs) %in% unique(has_dup_semf$semf))

  chosen_emfs = purrr::map_at(chosen_emfs, has_dup_emf, function(in_semf){
    peak_2_gemf = purrr::map_df(seq(1, nrow(in_semf)), function(in_row){
      data.frame(Sample_Peak = in_semf[in_row, "Sample_Peak"][[1]],
                 grouped_emf = in_semf[in_row, "grouped_emf"],
                 stringsAsFactors = FALSE)
    })
    bad_gemfs = dplyr::filter(peak_2_gemf, Sample_Peak %in% dup_peaks)
    good_gemfs = dplyr::filter(peak_2_gemf, !(grouped_emf %in% bad_gemfs$grouped_emf))

    if (nrow(good_gemfs) > 0) {
      return(choose_emf(all_gemfs[unique(good_gemfs$grouped_emf)], peak_mz, keep_ratio))
    } else {
      return(NULL)
    }

  })

  chosen_null = purrr::map_lgl(chosen_emfs, is.null)
  chosen_emfs = chosen_emfs[!chosen_null]
  names(chosen_emfs) = paste0("SEMF.", seq(1, length(chosen_emfs)))
  chosen_emfs
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
#' @param peak_mz M/Z information for each and every peak
#' @param keep_ratio the ratio used to keep other highly voted EMFs
#'
#' @return list
#' @export
merge_duplicate_semfs = function(chosen_emfs, all_gemfs, peak_mz, keep_ratio = 0.9){
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
  chosen_duplicates = purrr::map(merged_gemfs, function(.x){
    choose_emf(all_gemfs[unique(.x$grouped_emf)], peak_mz, keep_ratio)
  })

  all_chosen_emfs = c(chosen_nondup, chosen_duplicates)
  names(all_chosen_emfs) = paste0("SEMF.", seq(1, length(all_chosen_emfs)))

  remove_duplicates_across_semfs(all_chosen_emfs, all_gemfs, peak_mz, keep_ratio)

}

extract_emfs = function(chosen_emfs){
  all_samples = unique(unlist(purrr::map(chosen_emfs, ~ .x$Sample)))

  null_sample_matrix = matrix(as.character(NA), nrow = 1, ncol = length(all_samples))
  colnames(null_sample_matrix) = all_samples

  all_emfs = internal_map$map_function(chosen_emfs, function(use_emf){
    #message(i_emf)
    #use_emf = chosen_emfs[[i_emf]]
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

      list(peak_matrix = emf_matrix[nap_order, ], peak_info = emf_nap[nap_order, ])
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

#' extract IMFs from EMFs
#'
#' There are use cases where we don't actually want to do statistical testing at the EMF level, but rather at the level
#' of individual IMFs. In that case, we want to end up with a matrix of heights and M/Z for each IMF. See `Details` for
#' how we handle multiple EMFs.
#'
#' @param emfs list of EMFs
#' @param emf_info the info object for the list of EMFs
#'
#' @details For an EMF, the `peak_matrix` matrices are compared row by row, and for those that are identical, only one
#'   copy of the height and M/Z matrices will be returned. In the cases where there is disagreement, multiple lines
#'   will be returned.
#'
#' @return list
#' @export
#'
extract_imfs = function(emfs){
  all_compared = purrr::map(emfs, compare_emfs)
}

compare_emfs = function(emf){
  peak_matrices = purrr::map(emf, ~ .x$peak_matrix)

  # there is an odd property of some of these EMFs. In some samples, one EMF will get assignments, but
  # the other EMF will not, leading to differences in the peak matrices. But if both EMFs have the same
  # number of isotopologues, then they should have both been seen in the same samples. This code
  # is an attempt to rectify that situation, by checking that either the same peaks were assigned the IMF
  # in some samples, or that no peaks were observed. If there is a sample with discrepant peaks observed
  # across the IMFs, then this will break and fall back to the original peak matrices
  n_row = purrr::map_int(peak_matrices, ~ nrow(.x))
  max_row = max(n_row)
  peak_matrices2 = peak_matrices

  for (irow in seq_len(max_row)) {
    use_matrices = n_row <= irow
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
        warning("there was a mismatch!")
        return("FALSE")
      }
    })

    peak_matrices2 = purrr::map(peak_matrices2[use_matrices], function(.x){
      .x[irow, ] = match_or_na
      .x
    })

  }

  all_peaks = do.call(rbind, peak_matrices2)

  if ("FALSE" %in% all_peaks) {
    all_peaks = do.call(rbind, peak_matrices)
  }

  n_rows = purrr::map_int(peak_matrices, ~ nrow(.x))
  max_rows = max(n_rows)

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

  df_groups = data.frame(groups = group_imfs, row_index = seq(1, length(group_imfs)))
  split_groups = split(df_groups, df_groups$groups)

  single_index = purrr::map_df(split_groups, ~ .x[1, ])

  out_height = compare_heights[single_index$row_index, ]
  out_mz = compare_mz[single_index$row_index, ]

  out_info = purrr::map(split_groups, ~ compare_info[.x$row_index, ])

  list(height = out_height, mz = out_mz, info = out_info)
}
