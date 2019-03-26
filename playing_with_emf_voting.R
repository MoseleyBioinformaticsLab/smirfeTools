# helping me figure out labeled EMF voting
library(dplyr)
assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")

# I should extract all of the e-values, and see if I should be filtering out the 0.99, and see how the
# multiple mapping drops.
# I should also be asking how many of the grouped EMFs have different numbers of peaks across samples...

get_sample_emfs = function(sample_assignments, sample_id){
  sample_assignments = dplyr::mutate(sample_assignments, emf.adduct = paste0(complete_EMF, ".", adduct_IMF))

  peaks_by_emf = split(sample_assignments$PeakID, sample_assignments$emf.adduct) %>%
    purrr::map2_dfr(., names(.), function(.x, .y){
      peaks = unique(.x)
      tmp_frame = data.frame(emf = .y, peak_chr = paste(peaks, collapse = ","), stringsAsFactors = FALSE)
      tmp_frame$peaks = list(peaks)
      tmp_frame
    })

  grouped_emf = split(peaks_by_emf, peaks_by_emf$peak_chr)

  grouped_emf_peaks = purrr::map(grouped_emf, function(x){
    use_peaks = x$peaks[1][[1]]
    peak_info = purrr::map(x$emf, ~ dplyr::filter(sample_assignments, PeakID %in% use_peaks, emf.adduct %in% .x))
    list(emf = x$emf,
         peak_info = peak_info,
         peaks = use_peaks,
         peaks_chr = paste0(peak_info[[1]]$Sample[1], ".", use_peaks)
    )
  })

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
# work with a single sudo emf
single_sudo = sudo_emf_list[[2]]

grouped_emfs = all_gemfs[single_sudo$grouped_emf]

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
choose_emf = function(grouped_emfs, single_sudo){
  n_emfs = rle(sort(single_sudo$emf.adduct))
  # everybody only has one result, awesome!
  prop_emfs = n_emfs$lengths / sum(n_emfs$lengths)
  if (max(prop_emfs) > 0.8) {
    use_emf = n_emfs$values[(prop_emfs == max(prop_emfs))]
  } else {
    use_emf = n_emfs$values
  }
  if (length(use_emf) > 1) {
    # now see if we can resolve multiple EMFs based on their mass error and e-values
    emf_data = purrr::map2_df(grouped_emfs, names(grouped_emfs), function(.x, .y){
      purrr::map2_df(.x$peak_info, .y, extract_e_value_mass_error, best_evalue = TRUE)
    })
  }

}

observed_mz = purrr::map_df(assigned_data, ~ .x$data %>% dplyr::filter(Measurement %in% "ObservedMZ"))
all_evalues = purrr::map2_df(all_gemfs, names(all_gemfs), function(.x, .y){
  purrr::map2_df(.x$peak_info, .y, extract_e_value_mass_error, best_evalue = TRUE)
})

sudo_evalues = furrr::future_map_dfr(sudo_emf_list, function(.x){
  dplyr::filter(all_evalues, grouped_emf %in% .x$grouped_emf) %>%
    split(., .$grouped_emf) %>% purrr::map_df(., ~ dplyr::slice(.x, which.min(e_value)))
})

sudo_evalues_mz = dplyr::left_join(sudo_evalues, observed_mz[, c("Sample_Peak", "Value")], by = "Sample_Peak") %>%
  dplyr::mutate(e_value = as.numeric(e_value))

ggplot(sudo_evalues_mz, aes(x = Value, y = e_value)) + geom_point()
