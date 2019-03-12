# helping me figure out labeled EMF voting
library(dplyr)
assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")

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

  names(grouped_emf_peaks) = paste0("GEMF.", sample_id, ".", seq_along(grouped_emf_peaks))

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

# work with a single sudo emf
single_sudo = sudo_emf_list[[1]]
