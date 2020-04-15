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
json_files = gsub(".zip$", ".json", zip_files)
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
dplyr::group_by(coefficient_df, cluster) %>% dplyr::summarize(n = dplyr::n())

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

# Looking at low m/z frequency SDs using all the samples by cluster -----
# It turns out, that as soon as we have EMFs *across* clusters, the SD's in
# frequency space jump up drastically. It also turns out, that you can't move
# the peaks in frequency space at all by recomputing their location in frequency
# using a single model of M/Z to frequency. All of the points move by 1e-4 or 1e-5,
# which does not fix SDs on the order of 10 to 100.
#
# So, next possible steps:
# 1 - Combine within cluster SDs to find the actual SD, and then convert to M/Z
# 2 - We know that frequency can model the offset in M/Z quite well (need to check error),
#   can we use the relative change in M/Z by frequency, and the SD vs 0.5 offset ratio
#   to inflate the M/Z cutoff slightly?
# 3 - Worst case, do frequency SD in one group, work out the *actual* cutoff in M/Z,
#   and then use this M/Z cutoff across groups?
devtools::load_all("~/Projects/work/smirfeTools")
library(ggplot2)
library(dplyr)
assigned_data = readRDS("testing_stuff/lung_assigned_all.rds")
all_json = readRDS("testing_stuff/lung_json_all.rds")
coefficient_df = readRDS("testing_stuff/lung_sqrt_coefficients.rds")

coefficient_df = dplyr::mutate(coefficient_df, cluster =
                                 dplyr::case_when(
                                   sqrt < 29800000 ~ 1,
                                   dplyr::between(sqrt, 29802000, 29802600) ~ 2,
                                   sqrt > 29802600 ~ 3
                                 ))

all_freq_models = purrr::map(all_json, ~ .x$peak$frequency_mz$frequency_coefficients)
use_freq_model = all_freq_models[[1]]
freq_description = all_json[[1]]$peak$frequency_mz$frequency_fit_description

assigned_data2 = assigned_data
is_equal = purrr::map(assigned_data, function(.x){
  purrr::map(names(.x$scan_level$ObservedMZ), function(in_peak){
    old_freq = .x$scan_level$ObservedFrequency[[in_peak]]
    new_freq = FTMS.peakCharacterization:::predict_exponentials(.x$scan_level$ObservedMZ[[in_peak]], use_freq_model, freq_description)
    all.equal(old_freq, new_freq)
  })
})


low_evalue_cutoff = 0.1
low_mz_cutoff = 500
remove_elements = "S"

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

semf_ncluster = purrr::map_df(confident_sudo_emfs, function(.x){
  in_samples = strsplit(.x$grouped_EMF, ".", fixed = TRUE) %>% purrr::map_chr(., ~ .x[2]) %>% unique(.)
  tmp_cluster = dplyr::filter(coefficient_df, sample %in% in_samples)
  n_cluster = length(unique(tmp_cluster$cluster))
  data.frame(semf = .x$sudo_EMF[1], n_sample = length(in_samples), n_cluster = n_cluster, stringsAsFactors = FALSE)
})

sd_df = dplyr::left_join(sd_df, semf_ncluster, by = "semf")

ggplot(sd_df, aes(x = sd)) + geom_histogram(bins = 100) + facet_wrap(~ n_cluster, ncol = 1)
