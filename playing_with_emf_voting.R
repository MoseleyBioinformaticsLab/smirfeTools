# helping me figure out labeled EMF voting
library(dplyr)
library(ggplot2)
#assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")
assigned_data = readRDS("lung_matched_tissue_raw_smirfe_assignments_2018-09-07_imf_adduct_fusion1_cancer.rds")

# I should extract all of the e-values, and see if I should be filtering out the 0.99, and see how the
# multiple mapping drops.
# I should also be asking how many of the grouped EMFs have different numbers of peaks across samples...

get_sample_emfs = function(sample_assignments, sample_id, evalue_cutoff = 0.98){
  sample_assignments = dplyr::mutate(sample_assignments, emf.adduct = paste0(complete_EMF, ".", adduct_IMF))
  e_values = dplyr::filter(sample_assignments, Type %in% "e_value") %>% dplyr::mutate(e_value = as.numeric(Assignment_Data))
  clique_size = dplyr::filter(sample_assignments, Type %in% "clique_size") %>% dplyr::mutate(clique_size = as.integer(Assignment_Data))

  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$emf.adduct) %>%
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
    peak_info = purrr::map(x$emf, ~ dplyr::filter(sample_assignments, PeakID %in% use_peaks, emf.adduct %in% .x))
    grouped_evalues = purrr::map_dbl(x$emf, ~ dplyr::filter(e_values, PeakID %in% use_peaks, emf.adduct %in% .x) %>% dplyr::slice(which.min(e_value)) %>% dplyr::pull(e_value))
    clique_size = purrr::map_int(x$emf, ~ dplyr::filter(clique_size, PeakID %in% use_peaks, emf.adduct %in% .x) %>% dplyr::slice(1) %>% dplyr::pull(clique_size))

    keep_emf = grouped_evalues <= evalue_cutoff
    x = x[keep_emf, ]
    if (sum(keep_emf) > 0) {
      return(list(emf = x$emf,
           peak_info = peak_info[keep_emf],
           e_values = dplyr::mutate(x, e_value = grouped_evalues[keep_emf]),
           min_e_value = min(grouped_evalues[keep_emf]),
           clique_size = clique_size[1],
           peaks = use_peaks,
           peaks_chr = paste0(peak_info[[1]]$Sample[1], ".", use_peaks)
      ))
    } else {
      return(NULL)
    }

  })

  gep_null = purrr::map_lgl(grouped_emf_peaks, is.null)
  grouped_emf_peaks = grouped_emf_peaks[!gep_null]
  names(grouped_emf_peaks) = paste0("GEMF_", seq_along(grouped_emf_peaks), ".", sample_id)

  emf_by_emf = purrr::map_dfr(grouped_emf_peaks, function(x){
    purrr::map_dfr(x$peak_info, ~ data.frame(emf.adduct = .x$emf.adduct[1],
                                             isotopologue_EMF = unique(dplyr::filter(.x, Type %in% "isotopologue_EMF") %>% dplyr::pull(Assignment_Data)), stringsAsFactors = FALSE))
  }) %>% split(., .$isotopologue_EMF)

  multi_evidence_emf = purrr::map(emf_by_emf, function(x){
    out_data = NULL
    if (nrow(x) > 1) {
      adducts = strsplit(x$emf.adduct, ".", fixed = TRUE) %>% purrr::map_chr(~ .x[2])

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

library(furrr)
plan(multiprocess)
within_sample_emfs = furrr::future_map(assigned_data, function(.x){
  tmp_assign = dplyr::filter(.x$assignments, !grepl("S", complete_EMF))
  get_sample_emfs(tmp_assign, .x$sample)
})


all_gemf_emf_mapping = furrr::future_map_dfr(within_sample_emfs, function(x){
  purrr::map2_dfr(x$grouped_emf, names(x$grouped_emf), function(.x, .y){
    data.frame(grouped_emf = .y, emf.adduct = .x$emf, stringsAsFactors = FALSE)
  })
})

sudo_emf_list = vector("list", length(unique(all_gemf_emf_mapping$grouped_emf)))

i_sudo = 1
while (nrow(all_gemf_emf_mapping) > 0) {
  #message(i_sudo)
  tmp_emf_mapping = dplyr::filter(all_gemf_emf_mapping, emf.adduct %in% all_gemf_emf_mapping$emf.adduct[1])

  emf_iter = dplyr::left_join(tmp_emf_mapping, all_gemf_emf_mapping, by = "grouped_emf") %>%
    dplyr::transmute(grouped_emf = grouped_emf, emf.adduct = emf.adduct.y) %>% unique()

  adduct_iter = dplyr::left_join(emf_iter, all_gemf_emf_mapping, by = "emf.adduct") %>%
    dplyr::transmute(grouped_emf = grouped_emf.y, emf.adduct = emf.adduct) %>% unique()

  while (nrow(adduct_iter) != nrow(emf_iter)) {
    emf_iter = dplyr::left_join(adduct_iter, all_gemf_emf_mapping, by = "grouped_emf") %>%
      dplyr::transmute(grouped_emf = grouped_emf, emf.adduct = emf.adduct.y) %>% unique()

    adduct_iter = dplyr::left_join(emf_iter, all_gemf_emf_mapping, by = "emf.adduct") %>%
      dplyr::transmute(grouped_emf = grouped_emf.y, emf.adduct = emf.adduct) %>% unique()
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


emf_peak_mapping = purrr::map2_dfr(all_gemfs, names(all_gemfs), function(.x, .y){
  data.frame(emf = .y, peaks = .x$peaks_chr, stringsAsFactors = FALSE)
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


extract_e_value_mass_error = function(peak_info, gemf_id, best_evalue = FALSE){
  evalue_masserror = dplyr::filter(peak_info, Type %in% c("e_value", "mass_error", "PeakID", "Sample_Peak", "clique_size")) %>%
    tidyr::spread(Type, Assignment_Data) %>% dplyr::mutate(grouped_emf = gemf_id)

  if (best_evalue) {
    evalue_masserror = dplyr::slice(evalue_masserror, which.min(e_value))
  }

  evalue_masserror
}

# this function does not account for the fact that a peak may be mapped to several EMFs
# This is merely to choose amongst a number of EMFs that may be possible for a single set
# of peaks
choose_emf = function(grouped_emfs, keep_ratio = 0.9){
  grouped_evalues = purrr::map_df(grouped_emfs, ~ .x$e_values) %>% dplyr::mutate(information = 1 - e_value)

  # I think the way to incorporate multi-adduct evidence is to add some small constant value to each
  # of the samples where that EMF also had corresponding evidence. Either that, or do an initial vote
  # with loose ratio cutoff (say 0.5), then filter for corresponding in each sample, and do the vote
  # again, with something added to boost the corresponding.

  emf_votes = dplyr::group_by(grouped_evalues, emf) %>% dplyr::summarise(sum_information = sum(information)) %>%
    dplyr::mutate(max_ratio = sum_information / max(sum_information))
  keep_emf = dplyr::filter(emf_votes, max_ratio >= keep_ratio)

  gemf_emf = purrr::map2_dfr(grouped_emfs, names(grouped_emfs), function(.x, .y){
    tmp_df = data.frame(grouped_emf = .y, emf = .x$emf, stringsAsFactors = FALSE)
    tmp_df$peaks = list(.x$peaks)
    tmp_df
  })

  matching_keep_emf = dplyr::filter(gemf_emf, emf %in% keep_emf$emf)

}