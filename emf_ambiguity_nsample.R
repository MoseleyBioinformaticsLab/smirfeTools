# helping me figure out labeled EMF voting
library(smirfeTools)
library(furrr)
plan(multiprocess)
set_internal_map(furrr::future_map)
#assigned_data = readRDS("unlabeled_example_zaytseva_pdx.rds")
assigned_data = readRDS("lung_matched_tissue_raw_smirfe_assignments_2018-09-07_imf_adduct_fusion1_cancer.rds")
set.seed(20190418)
reps = seq(1, 10)
sample_sizes = round(seq(0.1, 0.9, 0.1) * length(assigned_data))

get_stats = function(extracted_data){
  emf_data = extracted_data$emfs

  emf_counts = purrr::map_df(emf_data, function(g_emf){
    n_emf = length(g_emf)
    n_sample = max(purrr::map_int(g_emf, function(i_emf){
      sum(colSums(!is.na(i_emf$height)) > 0)
    }))
    data.frame(n_emf = n_emf, n_sample = n_sample, n_used = ncol(i_emf$peak_matrix))
  })
}

full_set = extract_assigned_data(assigned_data) %>% get_stats()
full_set$i_rep = 1

mini_sets = purrr::cross2(reps, sample_sizes) %>%
  purrr::map_df(function(.x){
    use_samples = sample(length(assigned_data), .x[[2]])
    tmp = extract_assigned_data(assigned_data[use_samples]) %>% get_stats()
    tmp$i_rep = .x[[1]]
    tmp
  })

all_sets_emf_stats = rbind(full_set, mini_sets)

saveRDS(all_sets_emf_stats, file = "emf_ambiguity_nsample.rds")