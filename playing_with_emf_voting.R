# helping me figure out labeled EMF voting
library(dplyr)
library(ggplot2)
library(furrr)
plan(multiprocess)
#assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")
assigned_data = readRDS("lung_matched_tissue_raw_smirfe_assignments_2018-09-07_imf_adduct_fusion1_cancer.rds")

# I should extract all of the e-values, and see if I should be filtering out the 0.99, and see how the
# multiple mapping drops.
# I should also be asking how many of the grouped EMFs have different numbers of peaks across samples...

peak_mz = purrr::map_df(assigned_data, ~ dplyr::filter(.x$data, Measurement %in% "ObservedMZ"))

# gettting sample level emfs -----

get_sample_emfs = function(sample_assignments, sample_id, evalue_cutoff = 0.98, use_corroborating = TRUE){
  sample_assignments = dplyr::filter(sample_assignments, !grepl("S", complete_EMF))
  sample_assignments = dplyr::mutate(sample_assignments, emf_Adduct = paste0(complete_EMF, ".", adduct_IMF))
  e_values = dplyr::filter(sample_assignments, Type %in% "e_value") %>% dplyr::mutate(e_value = as.numeric(Assignment_Data))
  clique_size = dplyr::filter(sample_assignments, Type %in% "clique_size") %>% dplyr::mutate(clique_size = as.integer(Assignment_Data))

  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$emf_Adduct) %>%
    purrr::map2_dfr(., names(.), function(.x, .y){
      peaks = unique(.x)
      tmp_frame = data.frame(emf_Adduct = .y, PeakID_chr = paste(peaks, collapse = ","), stringsAsFactors = FALSE)
      tmp_frame$PeakID = list(peaks)
      tmp_frame
    })

  grouped_emf = split(peaks_by_emf, peaks_by_emf$PeakID_chr)

  grouped_emf_peaks = purrr::map(seq(grouped_emf), function(ige){
    #message(ige)
    x = grouped_emf[[ige]]
    use_peaks = x$PeakID[1][[1]]
    peak_info = purrr::map(x$emf_Adduct, ~ dplyr::filter(sample_assignments, PeakID %in% use_peaks, emf_Adduct %in% .x))
    names(peak_info) = x$emf_Adduct
    grouped_evalues = purrr::map_dbl(x$emf_Adduct, ~ dplyr::filter(e_values, PeakID %in% use_peaks, emf_Adduct %in% .x) %>% dplyr::slice(which.min(e_value)) %>% dplyr::pull(e_value))
    grouped_clique_size = purrr::map_int(x$emf_Adduct, ~ dplyr::filter(clique_size, PeakID %in% use_peaks, emf_Adduct %in% .x) %>% dplyr::slice(1) %>% dplyr::pull(clique_size))

    keep_emf = grouped_evalues <= evalue_cutoff
    x = x[keep_emf, ]
    if (sum(keep_emf) > 0) {
      return(list(emf_Adduct = x$emf_Adduct,
           peak_info = peak_info[keep_emf],
           e_values = dplyr::mutate(x, e_value = grouped_evalues[keep_emf], Sample = sample_id, complete_EMF = stringr::str_split_fixed(emf_Adduct, "\\.", 2)[,1], type = "primary"),
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
    purrr::map_dfr(x$peak_info, ~ data.frame(emf_Adduct = .x$emf_Adduct[1],
                                             isotopologue_EMF = unique(dplyr::filter(.x, Type %in% "isotopologue_EMF") %>% dplyr::pull(Assignment_Data)), stringsAsFactors = FALSE))
  }) %>% split(., .$isotopologue_EMF)

  multi_evidence_emf = purrr::map_df(emf_by_emf, function(x){
    out_data = NULL
    if (nrow(x) > 1) {
      #emf_adducts = stringr::str_split_fixed(x$emf_Adduct, "\\.", 2)
      adducts = strsplit(x$emf_Adduct, ".", fixed = TRUE) %>% purrr::map_chr(~ .x[2])

      if (sum(adducts %in% c("1H1", "14N1,1H4")) != length(adducts)) {
        out_data = x
        out_data
      }
    }
    out_data
  })

  if (use_corroborating) {
    all_evalues = purrr::map_df(grouped_emf_peaks, ~ .x$e_values)

    gep_2_emf_adduct = purrr::map2_dfr(grouped_emf_peaks, names(grouped_emf_peaks),
                                       function(.x, .y){
                                         #message(.y)
                                         data.frame(gep = .y, emf_Adduct = .x$e_values$emf_Adduct,
                                                    stringsAsFactors = FALSE)
                                       })
    gep_2_emf_adduct = dplyr::filter(gep_2_emf_adduct, emf_Adduct %in% multi_evidence_emf$emf_Adduct)
    unique_gep = unique(gep_2_emf_adduct$gep)
    has_multi_evidence = which(names(grouped_emf_peaks) %in% unique_gep)

    grouped_emf_peaks = purrr::map_at(grouped_emf_peaks, has_multi_evidence, function(gep){
      curr_evalues = gep$e_values

      multi_evalues = purrr::map_df(seq_len(nrow(curr_evalues)), function(curr_row){
        base_data = curr_evalues[curr_row, ]

        use_iso_emf = dplyr::filter(multi_evidence_emf, emf_Adduct %in% base_data$emf_Adduct) %>%
          dplyr::pull(isotopologue_EMF)
        emf_matches = dplyr::filter(multi_evidence_emf, !(emf_Adduct %in% base_data$emf_Adduct), isotopologue_EMF %in% use_iso_emf) %>%
          dplyr::pull(emf_Adduct)

        if (length(emf_matches) > 0) {
          other_data = dplyr::filter(all_evalues, emf_Adduct %in% emf_matches) %>%
            dplyr::mutate(emf_Adduct = base_data$emf_Adduct, type = "corroborating")
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

within_sample_emfs = furrr::future_map(assigned_data, function(.x){
  tmp_assign = dplyr::filter(.x$assignments, !grepl("S", complete_EMF))
  get_sample_emfs(tmp_assign, .x$sample, evalue_cutoff = 0.5)
})


all_gemf_emf_mapping = furrr::future_map_dfr(within_sample_emfs, function(x){
  purrr::map2_dfr(x$grouped_emf, names(x$grouped_emf), function(.x, .y){
    data.frame(grouped_emf = .y, emf_Adduct = .x$emf_Adduct, stringsAsFactors = FALSE)
  })
})

# sudo emf list generation -----
create_sudo_emfs = function(gemf_2_emf){
  sudo_emf_list = vector("list", length(unique(gemf_2_emf$grouped_emf)))

  i_sudo = 1
  while (nrow(gemf_2_emf) > 0) {
    #message(i_sudo)
    tmp_emf_mapping = dplyr::filter(gemf_2_emf, emf_Adduct %in% gemf_2_emf$emf_Adduct[1])

    emf_iter = dplyr::left_join(tmp_emf_mapping, gemf_2_emf, by = "grouped_emf") %>%
      dplyr::transmute(grouped_emf = grouped_emf, emf_Adduct = emf_Adduct.y) %>% unique()

    adduct_iter = dplyr::left_join(emf_iter, gemf_2_emf, by = "emf_Adduct") %>%
      dplyr::transmute(grouped_emf = grouped_emf.y, emf_Adduct = emf_Adduct) %>% unique()

    while (nrow(adduct_iter) != nrow(emf_iter)) {
      emf_iter = dplyr::left_join(adduct_iter, gemf_2_emf, by = "grouped_emf") %>%
        dplyr::transmute(grouped_emf = grouped_emf, emf_Adduct = emf_Adduct.y) %>% unique()

      adduct_iter = dplyr::left_join(emf_iter, gemf_2_emf, by = "emf_Adduct") %>%
        dplyr::transmute(grouped_emf = grouped_emf.y, emf_Adduct = emf_Adduct) %>% unique()
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

sudo_emf_list = create_sudo_emfs(all_gemf_emf_mapping)
# next things:
all_gemfs = unlist(purrr::map(within_sample_emfs, "grouped_emf"), recursive = FALSE)
# all_multi_evidence = purrr::map_dfr(within_sample_emfs, ~ purrr::map_dfr(.x$multi_evidence, ~.x)) %>%
#   split(., .$isotopologue_EMF)
#
#
# extract_e_value_mass_error = function(peak_info, gemf_id, best_evalue = FALSE){
#   evalue_masserror = dplyr::filter(peak_info, Type %in% c("e_value", "mass_error", "PeakID", "Sample_Peak", "clique_size")) %>%
#     tidyr::spread(Type, Assignment_Data) %>% dplyr::mutate(grouped_emf = gemf_id)
#
#   if (best_evalue) {
#     evalue_masserror = dplyr::slice(evalue_masserror, which.min(e_value))
#   }
#
#   evalue_masserror
# }
#
# emf_peak_mapping = purrr::map2_dfr(all_gemfs, names(all_gemfs), function(.x, .y){
#   data.frame(emf = .y, peaks = .x$Sample_Peak, stringsAsFactors = FALSE)
# })
#
# # look at the variation of e-value with respect to m/z
# observed_mz = purrr::map_df(assigned_data, ~ .x$data %>% dplyr::filter(Measurement %in% "ObservedMZ"))
# all_evalues = furrr::future_map2_dfr(all_gemfs, names(all_gemfs), function(.x, .y){
#   purrr::map2_df(.x$peak_info, .y, extract_e_value_mass_error, best_evalue = TRUE)
# })
#
# sudo_evalues = furrr::future_map_dfr(sudo_emf_list, function(.x){
#   dplyr::filter(all_evalues, grouped_emf %in% unique(.x$grouped_emf)) %>%
#     split(., .$grouped_emf) %>% purrr::map_df(., ~ dplyr::slice(.x, which.min(e_value)))
# })
#
# sudo_evalues_mz = dplyr::left_join(sudo_evalues, observed_mz[, c("Sample_Peak", "Value")], by = "Sample_Peak") %>%
#   dplyr::mutate(e_value = as.numeric(e_value))
#
# ggplot(sudo_evalues_mz, aes(x = Value, y = e_value)) + geom_point(alpha = 0.5)
#
# library(IRanges)
#
# mz_ranges = IRanges::IRanges(start = seq(floor(min(sudo_evalues_mz$Value)), ceiling(max(sudo_evalues_mz$Value)), 1) * 2e5,
#                              width = 10*2e5)
# mcols(mz_ranges)$mz = seq(floor(min(sudo_evalues_mz$Value)), ceiling(max(sudo_evalues_mz$Value)), 1)
#
# sudo_evalues_ranges = IRanges::IRanges(start = round(sudo_evalues_mz$Value * 2e5), width = 1)
# low_evalue = sudo_evalues_mz$e_value <= 0.1
#
# low_counts = IRanges::countOverlaps(mz_ranges, sudo_evalues_ranges[low_evalue])
# high_counts = IRanges::countOverlaps(mz_ranges, sudo_evalues_ranges[!low_evalue])
#
# mcols(mz_ranges)$low = low_counts + 1
# mcols(mz_ranges)$high = high_counts + 1
#
# count_mz_data = as.data.frame(mcols(mz_ranges))
# count_mz_data = dplyr::mutate(count_mz_data, norm = high + low, high2 = high / norm, low2 = low / norm,
#                               ratio = log(high / low), ratio2 = log(high2 / low2))
# count_mz_data_long = tidyr::gather(count_mz_data, key = "low_high", value = "count", low, high, low2, high2)
# ggplot(dplyr::filter(count_mz_data_long, low_high %in% c("low2", "high2")), aes(x = mz, y = count)) + geom_line() + facet_wrap(~ low_high, ncol = 1)
#
# ggplot(count_mz_data, aes(x = mz, y = ratio)) + geom_line()
#
# # what about the clique size? How does it vary across the sudo EMFs??
# clique_info = purrr::map_df(sudo_emf_list, function(x){
#   use_emfs = all_gemfs[unique(x$grouped_emf)]
#   emf_clique_sizes = purrr::map_int(use_emfs, "clique_size")
#   data.frame(mean = mean(emf_clique_sizes), sd = sd(emf_clique_sizes))
# })
#
# ggplot(clique_info, aes(x = mean, y = sd / mean))  + geom_jitter(alpha = 0.5)
#
# single_sudo = sudo_emf_list[[3]]
#
# grouped_emfs = all_gemfs[unique(single_sudo$grouped_emf)]
#
# by_sample = data.frame(grouped_emf = unique(single_sudo$grouped_emf), stringsAsFactors = FALSE) %>%
#   dplyr::mutate(gemf = stringr::str_split_fixed(grouped_emf, "\\.", 2)[,1],
#                 sample = stringr::str_split_fixed(grouped_emf, "\\.", 2)[,2]) %>%
#   dplyr::arrange(sample)

# voting on the emfs -----

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
                             emf_Adduct = imf_2_peak$emf_Adduct[1],
                             seq = tmp_mean$seq[1],
                             stringsAsFactors = FALSE)
    } else {
      match_imf = data.frame(complete_IMF = as.character(NA),
                             PeakID = NA,
                             Sample = as.character(NA),
                             Sample_peak = as.character(NA),
                             emf_Adduct = as.character(NA),
                             seq = as.integer(NA),
                             stringsAsFactors = FALSE)
    }
    match_imf
  })

  match_mz = match_mz[!is.na(match_mz$complete_IMF), ]

  null_evalue = data.frame(emf_Adduct = imf_2_peak$emf_Adduct[1], e_value = NA, stringsAsFactors = FALSE)
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

# this function does not account for the fact that a peak may be mapped to several EMFs
# This is merely to choose amongst a number of EMFs that may be possible for a single set
# of peaks
choose_emf = function(grouped_emfs, peak_mz, keep_ratio = 0.9){
  grouped_evalues = purrr::map_df(grouped_emfs, ~ .x$e_values) %>% dplyr::mutate(information = 1 - e_value)

  # I think the way to incorporate multi-adduct evidence is to add some small constant value to each
  # of the samples where that EMF also had corresponding evidence. Either that, or do an initial vote
  # with loose ratio cutoff (say 0.5), then filter for corresponding in each sample, and do the vote
  # again, with something added to boost the corresponding.

  emf_votes = dplyr::group_by(grouped_evalues, emf_Adduct) %>% dplyr::summarise(sum_information = sum(information)) %>%
    dplyr::mutate(max_ratio = sum_information / max(sum_information))
  keep_emf = dplyr::filter(emf_votes, max_ratio >= keep_ratio)

  null_evalue = data.frame(emf_Adduct = NA, e_value = NA, stringsAsFactors = FALSE)
  gemf_emf = purrr::map2_dfr(grouped_emfs, names(grouped_emfs), function(.x, .y){
    sample_from_grouped_emf = stringr::str_split_fixed(.y, "\\.", 2)[2]
    e_values = dplyr::filter(.x$e_values, emf_Adduct %in% keep_emf$emf_Adduct, type %in% "primary") %>% dplyr::select(emf_Adduct, e_value)

    if (nrow(e_values) > 0) {
      e_values$info = vector(mode = "list", length = nrow(e_values))
      e_values$Sample_Peak = vector(mode = "list", length = nrow(e_values))
      use_info = .x$peak_info[e_values$emf_Adduct]
      e_values$info = purrr::map(use_info, function(in_info){
        in_info %>% dplyr::filter(Type %in% "NAP") %>% dplyr::mutate(NAP = as.numeric(Assignment_Data)) %>%
          dplyr::select(complete_IMF, PeakID, Sample, Sample_Peak, emf_Adduct, NAP) %>%
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

      imf_by_emf = split(has_imf, has_imf$emf_Adduct)

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


chosen_emfs = purrr::map2(sudo_emf_list, names(sudo_emf_list), function(.x, .y){
 #message(.y)
 choose_emf(all_gemfs[unique(.x$grouped_emf)], peak_mz)
})

# metrics of voted sudo emfs -----

calculate_sudo_nemf = function(voted_emf_list){
  unique_emfs = unique(voted_emf_list$emf_Adduct)
  emf_adduct = stringr::str_split_fixed(unique_emfs, "\\.", 2)
  split_emfs = split(emf_adduct[,2], emf_adduct[,1])

  n_adduct = purrr::map_dbl(split_emfs, function(multi_adducts){
    ammonia_proton_mixup = c("14N1,1H4", "1H1")
    n_adduct = 1
    if (all(ammonia_proton_mixup %in% multi_adducts)) {
      multi_adducts = multi_adducts[!(multi_adducts %in% ammonia_proton_mixup)]
      n_adduct = 1
    }
    if (length(multi_adducts) > 0) {
      n_adduct = n_adduct + length(multi_adducts)
    }
    n_adduct
  })

  n_emf = sum(n_adduct)

  n_sample = length(unique(voted_emf_list$Sample))
  n_match_mz = length(unique(voted_emf_list[is.na(voted_emf_list$e_value), "sample"]))

  data.frame(n_emf = n_emf, n_sample = n_sample, n_matched_mz = n_match_mz)
}

voted_metrics = purrr::map_df(chosen_emfs, calculate_sudo_nemf)


peak_2_voted_emf = purrr::map2_df(chosen_emfs, names(chosen_emfs), function(.x, .y){
  all_peaks = unique(unlist(purrr::map(.x$Sample_Peak, ~ .x), use.names = FALSE))
  data.frame(Sample_Peak = all_peaks, semf = .y, stringsAsFactors = FALSE)
})

# statistics ----

##  Just using a cutoff of 0.98, and no corroborating information:
##  2258 sudo EMFs
##  1580 of 2258 have less than 5 samples associated
##  495 EMFs have more than 10 samples (of 50)
##  376 EMFs have a single EMF
##  1347 EMFs have <= 5 voted EMFs
##  4246 of 41K peaks are shared across EMFs (10%)

##  Cutoff of 0.98, with corroborating information
##  2258 sudo EMFs
##  1583 have less than 5 samples associated
##  489 EMFs have more than 10 samples
##  328 EMFs have a single EMF
##  1712 EMFs have <= 5 voted EMFs
##  4116 of  41K peaks are shared across EMFs (10%)

##  Cutoff of 0.75, with corroborating information
##  1925 sudo EMFs
##  1348 have less than 5 samples associated
##  431 EMFs have more than 10 samples
##  240 EMFs have a single EMF
##  1391 have <= 5 voted EMFs
##  1500 of 37K peaks are shared across EMFs (4%)

##  Cutoff of 0.50, with corroborating information
##  1746 sudo EMFs
##  1220 have less than 5 samples associated
##  409 have >= 10 samples
##  205 have a single EMF
##  1217 have <= 5 voted EMFs
##  903 of 35K peaks are shared across EMFs (2%)


nrow(voted_metrics)
dplyr::filter(voted_metrics, n_sample <= 5) %>% dplyr::summarize(n = n())
dplyr::filter(voted_metrics, n_sample >= 10) %>% dplyr::summarize(n = n())
dplyr::filter(voted_metrics, n_emf == 1) %>% dplyr::summarise(n = n())
dplyr::filter(voted_metrics, n_emf <= 5) %>% dplyr::summarise(n = n())
sum(duplicated(peak_2_voted_emf$Sample_Peak)) / nrow(peak_2_voted_emf)
sum(duplicated(peak_2_voted_emf$Sample_Peak))
nrow(peak_2_voted_emf)
