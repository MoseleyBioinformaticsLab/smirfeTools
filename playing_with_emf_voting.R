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

get_sample_emfs = function(sample_assignments, sample_id, evalue_cutoff = 0.98){
  sample_assignments = dplyr::mutate(sample_assignments, emf_Adduct = paste0(complete_EMF, ".", adduct_IMF))
  e_values = dplyr::filter(sample_assignments, Type %in% "e_value") %>% dplyr::mutate(e_value = as.numeric(Assignment_Data))
  clique_size = dplyr::filter(sample_assignments, Type %in% "clique_size") %>% dplyr::mutate(clique_size = as.integer(Assignment_Data))

  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$emf_Adduct) %>%
    purrr::map2_dfr(., names(.), function(.x, .y){
      peaks = unique(.x)
      tmp_frame = data.frame(emf = .y, peak_chr = paste(peaks, collapse = ","), stringsAsFactors = FALSE)
      tmp_frame$peaks = list(peaks)
      tmp_frame
    })

  grouped_emf = split(peaks_by_emf, peaks_by_emf$peak_chr)

  grouped_emf_peaks = purrr::map(seq(grouped_emf), function(ige){
    #message(ige)
    x = grouped_emf[[ige]]
    use_peaks = x$peaks[1][[1]]
    peak_info = purrr::map(x$emf, ~ dplyr::filter(sample_assignments, PeakID %in% use_peaks, emf_Adduct %in% .x))
    names(peak_info) = x$emf
    grouped_evalues = purrr::map_dbl(x$emf, ~ dplyr::filter(e_values, PeakID %in% use_peaks, emf_Adduct %in% .x) %>% dplyr::slice(which.min(e_value)) %>% dplyr::pull(e_value))
    clique_size = purrr::map_int(x$emf, ~ dplyr::filter(clique_size, PeakID %in% use_peaks, emf_Adduct %in% .x) %>% dplyr::slice(1) %>% dplyr::pull(clique_size))

    keep_emf = grouped_evalues <= evalue_cutoff
    x = x[keep_emf, ]
    if (sum(keep_emf) > 0) {
      return(list(emf = x$emf,
           peak_info = peak_info[keep_emf],
           e_values = dplyr::mutate(x, e_value = grouped_evalues[keep_emf]),
           min_e_value = min(grouped_evalues[keep_emf]),
           clique_size = clique_size[1],
           peaks = use_peaks,
           Sample_Peak = paste0(peak_info[[1]]$Sample[1], "_", use_peaks)
      ))
    } else {
      return(NULL)
    }

  })

  gep_null = purrr::map_lgl(grouped_emf_peaks, is.null)
  grouped_emf_peaks = grouped_emf_peaks[!gep_null]
  names(grouped_emf_peaks) = paste0("GEMF_", seq_along(grouped_emf_peaks), ".", sample_id)

  emf_by_emf = purrr::map_dfr(grouped_emf_peaks, function(x){
    purrr::map_dfr(x$peak_info, ~ data.frame(emf_Adduct = .x$emf_Adduct[1],
                                             isotopologue_EMF = unique(dplyr::filter(.x, Type %in% "isotopologue_EMF") %>% dplyr::pull(Assignment_Data)), stringsAsFactors = FALSE))
  }) %>% split(., .$isotopologue_EMF)

  multi_evidence_emf = purrr::map(emf_by_emf, function(x){
    out_data = NULL
    if (nrow(x) > 1) {
      adducts = strsplit(x$emf_Adduct, ".", fixed = TRUE) %>% purrr::map_chr(~ .x[2])

      if (sum(adducts %in% c("1H1", "14N1,1H4")) != length(adducts)) {
        out_data = x
        out_data$sample_id = sample_id
        out_data
      }
    }
    out_data
  })

  multi_evidence_emf = multi_evidence_emf[purrr::map_lgl(multi_evidence_emf, ~ !is.null(.x))]

  return(list(grouped_emf = grouped_emf_peaks, multi_evidence = multi_evidence_emf))
  # we should add the peaks to GEMF mapping here, so we can do what is noted on the last line.
}

within_sample_emfs = furrr::future_map(assigned_data, function(.x){
  tmp_assign = dplyr::filter(.x$assignments, !grepl("S", complete_EMF))
  get_sample_emfs(tmp_assign, .x$sample)
})


all_gemf_emf_mapping = furrr::future_map_dfr(within_sample_emfs, function(x){
  purrr::map2_dfr(x$grouped_emf, names(x$grouped_emf), function(.x, .y){
    data.frame(grouped_emf = .y, emf_Adduct = .x$emf, stringsAsFactors = FALSE)
  })
})

sudo_emf_list = vector("list", length(unique(all_gemf_emf_mapping$grouped_emf)))

i_sudo = 1
while (nrow(all_gemf_emf_mapping) > 0) {
  #message(i_sudo)
  tmp_emf_mapping = dplyr::filter(all_gemf_emf_mapping, emf_Adduct %in% all_gemf_emf_mapping$emf_Adduct[1])

  emf_iter = dplyr::left_join(tmp_emf_mapping, all_gemf_emf_mapping, by = "grouped_emf") %>%
    dplyr::transmute(grouped_emf = grouped_emf, emf_Adduct = emf_Adduct.y) %>% unique()

  adduct_iter = dplyr::left_join(emf_iter, all_gemf_emf_mapping, by = "emf_Adduct") %>%
    dplyr::transmute(grouped_emf = grouped_emf.y, emf_Adduct = emf_Adduct) %>% unique()

  while (nrow(adduct_iter) != nrow(emf_iter)) {
    emf_iter = dplyr::left_join(adduct_iter, all_gemf_emf_mapping, by = "grouped_emf") %>%
      dplyr::transmute(grouped_emf = grouped_emf, emf_Adduct = emf_Adduct.y) %>% unique()

    adduct_iter = dplyr::left_join(emf_iter, all_gemf_emf_mapping, by = "emf_Adduct") %>%
      dplyr::transmute(grouped_emf = grouped_emf.y, emf_Adduct = emf_Adduct) %>% unique()
  }

  sudo_emf_list[[i_sudo]] = adduct_iter
  all_gemf_emf_mapping = dplyr::filter(all_gemf_emf_mapping, !(grouped_emf %in% adduct_iter$grouped_emf))
  i_sudo = i_sudo + 1
}

keep_sudo = purrr::map_lgl(sudo_emf_list, ~ !is.null(.x))
sudo_emf_list = sudo_emf_list[keep_sudo]

# next things:
all_gemfs = unlist(purrr::map(within_sample_emfs, "grouped_emf"), recursive = FALSE)
all_multi_evidence = purrr::map_dfr(within_sample_emfs, ~ purrr::map_dfr(.x$multi_evidence, ~.x)) %>%
  split(., .$isotopologue_EMF)


extract_e_value_mass_error = function(peak_info, gemf_id, best_evalue = FALSE){
  evalue_masserror = dplyr::filter(peak_info, Type %in% c("e_value", "mass_error", "PeakID", "Sample_Peak", "clique_size")) %>%
    tidyr::spread(Type, Assignment_Data) %>% dplyr::mutate(grouped_emf = gemf_id)

  if (best_evalue) {
    evalue_masserror = dplyr::slice(evalue_masserror, which.min(e_value))
  }

  evalue_masserror
}

emf_peak_mapping = purrr::map2_dfr(all_gemfs, names(all_gemfs), function(.x, .y){
  data.frame(emf = .y, peaks = .x$Sample_Peak, stringsAsFactors = FALSE)
})

# look at the variation of e-value with respect to m/z
observed_mz = purrr::map_df(assigned_data, ~ .x$data %>% dplyr::filter(Measurement %in% "ObservedMZ"))
all_evalues = furrr::future_map2_dfr(all_gemfs, names(all_gemfs), function(.x, .y){
  purrr::map2_df(.x$peak_info, .y, extract_e_value_mass_error, best_evalue = TRUE)
})

sudo_evalues = furrr::future_map_dfr(sudo_emf_list, function(.x){
  dplyr::filter(all_evalues, grouped_emf %in% unique(.x$grouped_emf)) %>%
    split(., .$grouped_emf) %>% purrr::map_df(., ~ dplyr::slice(.x, which.min(e_value)))
})

sudo_evalues_mz = dplyr::left_join(sudo_evalues, observed_mz[, c("Sample_Peak", "Value")], by = "Sample_Peak") %>%
  dplyr::mutate(e_value = as.numeric(e_value))

ggplot(sudo_evalues_mz, aes(x = Value, y = e_value)) + geom_point(alpha = 0.5)

library(IRanges)

mz_ranges = IRanges::IRanges(start = seq(floor(min(sudo_evalues_mz$Value)), ceiling(max(sudo_evalues_mz$Value)), 1) * 2e5,
                             width = 10*2e5)
mcols(mz_ranges)$mz = seq(floor(min(sudo_evalues_mz$Value)), ceiling(max(sudo_evalues_mz$Value)), 1)

sudo_evalues_ranges = IRanges::IRanges(start = round(sudo_evalues_mz$Value * 2e5), width = 1)
low_evalue = sudo_evalues_mz$e_value <= 0.1

low_counts = IRanges::countOverlaps(mz_ranges, sudo_evalues_ranges[low_evalue])
high_counts = IRanges::countOverlaps(mz_ranges, sudo_evalues_ranges[!low_evalue])

mcols(mz_ranges)$low = low_counts + 1
mcols(mz_ranges)$high = high_counts + 1

count_mz_data = as.data.frame(mcols(mz_ranges))
count_mz_data = dplyr::mutate(count_mz_data, norm = high + low, high2 = high / norm, low2 = low / norm,
                              ratio = log(high / low), ratio2 = log(high2 / low2))
count_mz_data_long = tidyr::gather(count_mz_data, key = "low_high", value = "count", low, high, low2, high2)
ggplot(dplyr::filter(count_mz_data_long, low_high %in% c("low2", "high2")), aes(x = mz, y = count)) + geom_line() + facet_wrap(~ low_high, ncol = 1)

ggplot(count_mz_data, aes(x = mz, y = ratio)) + geom_line()

# what about the clique size? How does it vary across the sudo EMFs??
clique_info = purrr::map_df(sudo_emf_list, function(x){
  use_emfs = all_gemfs[unique(x$grouped_emf)]
  emf_clique_sizes = purrr::map_int(use_emfs, "clique_size")
  data.frame(mean = mean(emf_clique_sizes), sd = sd(emf_clique_sizes))
})

ggplot(clique_info, aes(x = mean, y = sd / mean))  + geom_jitter(alpha = 0.5)

single_sudo = sudo_emf_list[[3]]

grouped_emfs = all_gemfs[unique(single_sudo$grouped_emf)]

by_sample = data.frame(grouped_emf = unique(single_sudo$grouped_emf), stringsAsFactors = FALSE) %>%
  dplyr::mutate(gemf = stringr::str_split_fixed(grouped_emf, "\\.", 2)[,1],
                sample = stringr::str_split_fixed(grouped_emf, "\\.", 2)[,2]) %>%
  dplyr::arrange(sample)

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

  unknown_mz = dplyr::left_join(unknown_peaks, peak_mz, by = c("Peaks" = "Sample_Peak"))
  # this also needs to have the peaks sorted by their NAP, and match in order
  # of NAP. When the next one doesn't match, stop!
  match_mz = purrr::map_df(mean_mz$complete_IMF, function(test_imf){
    message(test_imf)
    test_mz = dplyr::filter(mean_mz, complete_IMF %in% test_imf) %>%
      dplyr::slice(rep(1:n()), each = 3)
    test_mz = test_mz[rep(1, each = nrow(unknown_mz)), ]
    tmp_mean = test_mz %>% dplyr::mutate(diff = abs(mean - unknown_mz$Value),
                                         Sample_Peak = unknown_mz$Peaks,
                                         PeakID = unknown_mz$PeakID,
                                         Sample = unknown_mz$Sample.x) %>%
      dplyr::filter(diff <= sd2)

    if (nrow(tmp_mean) > 0) {
      tmp_mean = dplyr::slice(tmp_mean, which.min(diff))
      match_imf = data.frame(complete_IMF = tmp_mean$complete_IMF,
                             PeakID = tmp_mean$PeakID,
                             Sample = tmp_mean$Sample,
                             Sample_Peak = tmp_mean$Sample_Peak,
                             emf_Adduct = imf_2_peak$emf_Adduct[1],
                             seq = tmp_mean$seq,
                             stringsAsFactors = FALSE)

    } else {
      match_imf = data.frame(complete_IMF = as.character(NA),
                             PeakID = as.integer(NA),
                             Sample = as.character(NA),
                             Sample_Peak = as.character(NA),
                             emf_Adduct = as.character(NA),
                             seq = as.integer(NA),
                             stringsAsFactors = FALSE)
    }
    match_imf
  })

  match_mz = match_mz[!is.na(match_mz$complete_IMF), ]
  print(match_mz)
  null_evalue = data.frame(emf = imf_2_peak$emf_Adduct[1], e_value = NA, stringsAsFactors = FALSE)
  if (nrow(match_mz) > 0) {
    seq_match = all(match_mz$seq == seq_len(nrow(match_mz)))
    unique_match = length(unique(match_mz$complete_IMF)) == nrow(match_mz)

    if (!(seq_match && unique_match)) {
      match_mz = match_mz[1,]
      match_mz[1, ] = NA
    }
  } else {
    match_mz = match_mz[1,]
    match_mz[1, ] = NA
  }

  out_evalues = null_evalue
  out_evalues$info = list(match_mz)
  out_evalues$peaks = list(match_mz$Sample_Peak)
  out_evalues$sample = match_mz$Sample[1]
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

  emf_votes = dplyr::group_by(grouped_evalues, emf) %>% dplyr::summarise(sum_information = sum(information)) %>%
    dplyr::mutate(max_ratio = sum_information / max(sum_information))
  keep_emf = dplyr::filter(emf_votes, max_ratio >= keep_ratio)

  null_evalue = data.frame(emf = NA, e_value = NA, stringsAsFactors = FALSE)
  gemf_emf = purrr::map2_dfr(grouped_emfs, names(grouped_emfs), function(.x, .y){
    sample_from_grouped_emf = stringr::str_split_fixed(.y, "\\.", 2)[2]
    e_values = dplyr::filter(.x$e_values, emf %in% keep_emf$emf) %>% dplyr::select(emf, e_value)

    if (nrow(e_values) > 0) {
      e_values$info = vector(mode = "list", length = nrow(e_values))
      e_values$peaks = vector(mode = "list", length = nrow(e_values))
      use_info = .x$peak_info[e_values$emf]
      e_values$info = purrr::map(use_info, function(in_info){
        in_info %>% dplyr::filter(Type %in% "NAP") %>% dplyr::mutate(NAP = as.numeric(Assignment_Data)) %>%
          dplyr::select(complete_IMF, PeakID, Sample, Sample_Peak, emf_Adduct, NAP) %>%
          unique()
      })

      e_values$peaks = list(.x$Sample_Peak)
      e_values$sample = as.character(sample_from_grouped_emf)
      e_values$grouped_emf = as.character(.y)
    } else {
      e_values = null_evalue
      e_values$info = vector(mode = "list", length = nrow(e_values))
      e_values$peaks = list(.x$Sample_Peak)
      e_values$sample = as.character(sample_from_grouped_emf)
      e_values$grouped_emf = as.character(.y)
    }
    e_values
  })
  na_evalues = purrr::map_lgl(gemf_emf$e_value, is.na)

  have_imf_peaks = unique(unlist(purrr::map(gemf_emf$peaks[!na_evalues], ~ .x), use.names = FALSE))

  if (sum(na_evalues) > 0) {
    trimmed_gemf_emf = gemf_emf[!na_evalues, ]
    missing_imf = purrr::map(which(na_evalues), function(.x){
      tmp_df = data.frame(Peaks = as.character(gemf_emf$peaks[[.x]]), Sample = gemf_emf$sample[.x],
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

  if (length(setdiff(unique(gemf_emf$sample), unique(out_gemf_emf$sample))) > 0) {
    message("sample doesn't have IMFs")
  }

  out_gemf_emf

}

names(sudo_emf_list) = paste0("SEMF.", seq(1, length(sudo_emf_list)))
chosen_emfs = purrr::map2(sudo_emf_list[1:100], names(sudo_emf_list)[1:100], function(.x, .y){
 message(.y)
 choose_emf(all_gemfs[unique(.x$grouped_emf)], peak_mz)
})
#
# # debugging
grouped_emfs = all_gemfs[unique(sudo_emf_list[[16]]$grouped_emf)]

# chosen_emfs = furrr::future_map(sudo_emf_list, ~ choose_emf(all_gemfs[unique(.x$grouped_emf)]))
