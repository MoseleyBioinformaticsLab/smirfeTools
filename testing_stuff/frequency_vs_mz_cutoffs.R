library(smirfeTools)
library(knitrProgressBar)
library(furrr)
plan(multiprocess)
set_internal_map(furrr::future_map)
options(future.globals.maxSize = 800 * 1024 ^ 2)

save_dir = "/mlab/scratch/rmflight/smirfeTools_data"
setwd(save_dir)

assigned_files <- dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2020-01-27",
                      full.names = TRUE, pattern = "assigned.json")
pb <- progress_estimated(length(assigned_files))

assigned_data <- purrr::map(assigned_files, function(in_file){
  tmp_assign = read_smirfe_assignment(in_file)
  filter_data = score_filter_assignments(tmp_assign, filter_conditions = e_value <= 0.5)
  update_progress(pb)
  filter_data
})

names(assigned_data) = purrr::map_chr(assigned_data, ~ .x$sample)

all_emfs = unique(unlist(purrr::map(assigned_data, ~ .x$assignments$isotopologue_EMF)))


zip_files <- dir("/mlab/scratch/cesb_data/zip_files/lung_matched_tissue-2020-01-27",
                 full.names = TRUE, pattern = ".zip$")
all_coefficients = extract_coefficient_data(zip_files)

names(all_coefficients) = purrr::map_chr(all_coefficients, "sample")


coefficient_df = purrr::map_df(all_coefficients, function(in_sample){
  data.frame(sample = in_sample$sample,
             sqrt = in_sample$coefficients$frequency_coefficients[2],
             stringsAsFactors = FALSE)
})
saveRDS(coefficient_df, file = "lung_sqrt_coefficients.rds")
coefficient_cluster = kmeans(coefficient_df$sqrt, 3)
coefficient_df$cluster = coefficient_cluster$cluster

coefficients_by_cluster = split(coefficient_df, coefficient_df$cluster)

sd_by_cluster = purrr::map(coefficients_by_cluster, function(in_cluster){
  tmp_sd = find_confident_frequency_sd(assigned_data[in_cluster$sample])
  tmp_coefficients = all_coefficients[in_cluster$sample]

  use_sd = mean(tmp_sd$sd) * 2
  mz_diff = create_mz_diffs(tmp_coefficients, use_sd)
  list(mz = mz_diff, sd = tmp_sd, cluster = in_cluster$cluster[1])
})

mz_values = purrr::map_dbl(sd_by_cluster, ~ max(.x$mz$Value))
use_mz = sd_by_cluster[[which.max(mz_values)]]$mz

classified_emfs = import_emf_classifications("all_emfs_2020-01-27_classified.json")


lipid_df = weight_lipid_classifications(classified_emfs, lipid_weight = 2, not_lipid_weight = 2)
new_scores = purrr::map(assigned_data, function(in_data){
  score_filter_assignments(in_data, filter_conditions = e_value <= 0.5,
                           emf_weight = lipid_df)
})
new_extracts = extract_assigned_data(new_scores, difference_cutoff = use_mz,
                                     difference_measure = "ObservedMZ",
                                     progress = FALSE)

# try everything with cluster 1
saveRDS(sd_by_cluster, file = "lung_sd_by_cluster.rds")
saveRDS(coef)
cluster_id = 3
use_cluster = sd_by_cluster[[cluster_id]]
tmp_sd = use_cluster$sd
freq_sd = mean(tmp_sd$sd) * 2
new_scores_c1 = new_scores[coefficients_by_cluster[[cluster_id]]$sample]
new_extract_freq = extract_assigned_data(new_scores_c1, difference_cutoff = freq_sd,
                                        difference_measure = "ObservedFrequency",
                                        progress = TRUE)

tmp_coefficients = all_coefficients[coefficients_by_cluster[[cluster_id]]$sample]
mz_diff = create_mz_diffs(tmp_coefficients,freq_sd, mz_values = seq(150, 1600, 0.1))
new_extract_mz = extract_assigned_data(new_scores_c1, difference_cutoff = mz_diff,
                                       difference_measure = "ObservedMZ",
                                       progress = TRUE)

# this didn't work, we didn't get any trimming on the chosen_EMFs. So lets modify our function a bit
create_mz_diffs_local = function(coefficient_info, frequency_offset = NULL, mz_values = seq(150, 1600, 0.5)){
  if (is.null(frequency_offset)) {
    stop("frequency_offset is NULL, please supply a value!")
  }

  freq_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$frequency_coefficients
  }) %>% do.call(rbind, .)
  mz_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$mz_coefficients
  }) %>% do.call(rbind, .)

  # we are going to choose 1 specific model, in this case the one that has the
  # smallest sqrt model
  min_freq = which.max(freq_coefficients[, 2])
  use_freq = freq_coefficients[min_freq, ]
  use_mz = mz_coefficients[min_freq, ]

  freq_description = coefficient_info[[1]]$coefficients$frequency_fit_description
  mz_description = coefficient_info[[1]]$coefficients$mz_fit_description

  mz_2_freq = FTMS.peakCharacterization:::predict_exponentials(mz_values, use_freq, freq_description)

  frequency_diff = mz_2_freq + frequency_offset

  mz_init = FTMS.peakCharacterization:::predict_exponentials(mz_2_freq, use_mz, mz_description)
  mz_shift = FTMS.peakCharacterization:::predict_exponentials(frequency_diff, use_mz, mz_description)
  mz_diff = abs(mz_shift - mz_init)

  data.frame(Index = mz_values, Value = mz_diff)

}

mz_diff2 = create_mz_diffs_local(tmp_coefficients, freq_sd/11, mz_values = seq(150, 1600, 0.1))
new_extract_mz2 = extract_assigned_data(new_scores_c1, difference_cutoff = mz_diff2,
                                       difference_measure = "ObservedMZ",
                                       progress = TRUE)

new_extract_mz2$diff = mz_diff2
saveRDS(new_extract_mz2, file = "lung_mz2_based.rds")
saveRDS(new_extract_freq, file = "lung_freq_based.rds")


# ---- viz sd
sd_by_cluster = readRDS("testing_stuff/lung_sd_by_cluster.rds")


sd_df = purrr::map_df(sd_by_cluster, function(.x){
  tmp_sd = .x$sd
  tmp_sd$cluster = .x$cluster
  tmp_sd
})

library(ggplot2)
library(smirfeTools)
ggplot(sd_df, aes(x = sd, fill = as.factor(cluster))) + geom_histogram(bins = 100) + facet_wrap(~ cluster, ncol = 1)
ggplot(sd_df, aes(x = sd)) + geom_histogram(bins = 100) + geom_vline(xintercept = median(sd_df$sd), color = "red")
ggplot(sd_df, aes(x = sd)) + geom_density() + geom_vline(xintercept = median(sd_df$sd), color = "red")

coefficient_df = readRDS("testing_stuff/lung_sqrt_coefficients.rds")
ggplot(coefficient_df, aes(x = sqrt, fill = as.factor(cluster))) + geom_histogram(bins = 100)

lung_mz2_based = readRDS("testing_stuff/lung_mz2_based.rds")
lung_freq_based = readRDS("testing_stuff/lung_freq_based.rds")
mz_info = extract_imf_emf_data(lung_mz2_based$emfs, by = "IMF")

freq_info = extract_imf_emf_data(lung_freq_based$emfs, by = "IMF")
freq_loc = freq_info$location
freq_greater3 = rowSums(!is.na(freq_loc)) >= 3
freq_loc = freq_loc[freq_greater3, ]
freq_sd = data.frame(Index = apply(freq_loc, 1, mean, na.rm = TRUE),
                     Value = apply(freq_loc, 1, sd, na.rm = TRUE),
                     type = "observed")
mz_diff = lung_mz2_based$diff
mz_diff$type = "predicted"
all_diff = rbind(freq_sd, mz_diff)

theme_set(cowplot::theme_cowplot())
ggplot(all_diff, aes(x = Index, y = Value, color = type)) + geom_point() +
  theme(legend.position = c(0.1, 0.7)) +
  labs(x = "M/Z", subtitle = "Comparing the M/Z based cutoff divided by 11 to the actual \nSD observed from using Frequency based cutoff")

mz_loc = mz_info$location
mz_greater3 = rowSums(!is.na(mz_loc)) >= 3
mz_loc = mz_loc[mz_greater3, ]
mz_sd = data.frame(Index = apply(mz_loc, 1, mean, na.rm = TRUE),
                   Value = apply(mz_loc, 1, sd, na.rm = TRUE),
                   type = "observed")
all_mz = rbind(mz_sd, mz_diff)
ggplot(all_mz, aes(x = Index, y = Value, color = type)) + geom_point() +
  theme(legend.position = c(0.1, 0.7)) +
  labs(x = "M/Z", subtitle = "Comparing the M/Z based cutoff divided by 11 to the actual\n M/Z SD observed from using M/Z based cutoff")

mz_org = readRDS("testing_stuff/lung_mzdiff_freqsd.rds")
mz_org$type = "original"

all_freq = rbind(mz_org, freq_sd)
ggplot(all_freq, aes(x = Index, y = Value, color = type)) + geom_point() +
  theme(legend.position = c(0.1, 0.7)) +
  labs(x = "M/Z", subtitle = "Comparing the directly predicted M/Z based cutoff to the actual\n M/Z SD observed when using the Frequency cutoff")


create_mz_diffs_frequency = function(coefficient_info, frequency_offset = NULL, assigned_info){
  if (is.null(frequency_offset)) {
    stop("frequency_offset is NULL, please supply a value!")
  }

  freq_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$frequency_coefficients
  }) %>% do.call(rbind, .)
  mz_coefficients = purrr::map(coefficient_info, function(in_sample){
    in_sample$coefficients$mz_coefficients
  }) %>% do.call(rbind, .)

  # we are going to choose 1 specific model, in this case the one that has the
  # smallest sqrt model
  min_freq = which.max(freq_coefficients[, 2])
  use_freq = freq_coefficients[min_freq, ]
  use_mz = mz_coefficients[min_freq, ]

  freq_description = coefficient_info[[1]]$coefficients$frequency_fit_description
  mz_description = coefficient_info[[1]]$coefficients$mz_fit_description

  use_sample = rownames(freq_coefficients)[min_freq]
  frequency_range = range(assigned_info[[use_sample]]$data$ObservedFrequency)
  frequency_points = seq(round(frequency_range[1]), round(frequency_range[2]), 1)
  frequency_off = frequency_points + frequency_offset

  mz_init = FTMS.peakCharacterization:::predict_exponentials(frequency_points, use_mz, mz_description)
  mz_shift = FTMS.peakCharacterization:::predict_exponentials(frequency_off, use_mz, mz_description)

  mz_diff = abs(mz_shift - mz_init)

  data.frame(Index = mz_init, Value = mz_diff)

}

assigned_data = readRDS("testing_stuff/lung_assigned_data.rds")
load("testing_stuff/lung_coefficient_info.rda")
frequency_offset = 1.160158
cluster_id = 3
use_coefficients = dplyr::filter(coefficient_df, cluster == cluster_id)
coefficient_info = all_coefficients[use_coefficients$sample]
assigned_info = assigned_data[use_coefficients$sample]

new_diff = create_mz_diffs_frequency(coefficient_info, freq_sd, cluster_assigned)
