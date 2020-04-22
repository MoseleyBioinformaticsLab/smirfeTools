# Example Workflow Across 3 Groups
library(smirfeTools)
library(ggplot2)
library(dplyr)
library(furrr)
library(metabolomicsUtilities)
theme_set(cowplot::theme_cowplot())
plan(multiprocess)
set_internal_map(furrr::future_map)
options(future.globals.maxSize = 800 * 1024 ^ 2)
assigned_data = readRDS("testing_stuff/lung_assigned_all.rds")
all_json = readRDS("testing_stuff/lung_json_all.rds")
coefficient_df = readRDS("testing_stuff/lung_sqrt_coefficients.rds")


coefficient_df = dplyr::mutate(coefficient_df, cluster =
                                 dplyr::case_when(
                                   sqrt < 29800000 ~ 1,
                                   dplyr::between(sqrt, 29802000, 29802600) ~ 2,
                                   sqrt > 29802600 ~ 3
                                 ))

coefficients_by_cluster = split(coefficient_df, coefficient_df$cluster)

sd_by_cluster = purrr::map_df(coefficients_by_cluster, function(in_cluster){
  tmp_sd = find_confident_frequency_sd(assigned_data[in_cluster$sample])
  tmp_sd$cluster = in_cluster$cluster[1]
  tmp_sd
})

mode_function = function(values){
  density_estimate = stats::density(values)
  mode_value = density_estimate$x[which.max(density_estimate$y)]
  mode_value
}

sd_mode = mode_function(sd_by_cluster$sd) * 2

classified_emfs = import_emf_classifications("testing_stuff/all_emfs_2020-01-27_classified.json")
lipid_df = weight_lipid_classifications(classified_emfs, lipid_weight = 2, not_lipid_weight = 2)
assigned_data = purrr::map(assigned_data, function(in_data){
  score_filter_assignments(in_data, filter_conditions = e_value <= 0.5,
                           emf_weight = lipid_df)
})
assigned_data = assigned_data[coefficient_df$sample]
split_assigned = split(assigned_data, coefficient_df$cluster)

get_mz_after_freq = function(assigned_set, freq_cutoff, set_id = NULL){
  extracts = extract_assigned_data(assigned_set, difference_cutoff = freq_cutoff,
                                   difference_measure = "ObservedFrequency",
                                   progress = FALSE)

  imf_mz = extract_imf_emf_data(extracts, Height, ObservedMZ, by = "IMF", scanlevel = TRUE)

  scan_mz = purrr::map_df(imf_mz$scanlevel_location, function(.x){
    data.frame(Index = mean(.x, na.rm = TRUE), Value = sd(.x, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })

  if (!is.null(set_id)) {
    scan_mz$set = set_id
  }
  scan_mz
}

mz_after_freq = purrr::imap_dfr(split_assigned, ~ get_mz_after_freq(.x, sd_mode, .y))

ggplot(dplyr::filter(mz_after_freq, Value <= 0.1), aes(x = Index, y = Value)) + geom_point() + facet_wrap(~ set, scales = "free", ncol = 1)

ggplot(dplyr::filter(mz_after_freq, Value <= 0.1), aes(x = Index, y = Value, color = as.factor(set))) + geom_point(alpha = 0.2)

mz_after_freq = dplyr::filter(mz_after_freq, Value <= 0.1)

mz_fit = mgcv::gam(Value ~ s(Index, bs = "cs"), data = mz_after_freq)
mz_after_freq$fit = mz_fit$fitted.values
mz_after_freq$fit2 = 2 * mz_after_freq$fit

ggplot(mz_after_freq, aes(x = Index, y = Value, color = as.factor(set))) + geom_point(alpha = 0.2) +
  geom_line(data = mz_after_freq, aes(x = Index, y = fit), color = "red") +
  geom_line(data = mz_after_freq, aes(x = Index, y = fit2), color = "blue")

# cool, looks like we can combine the fit across the sets
mz_index = data.frame(Index = seq(140, 1610, 0.5))
mz_pred = predict(mz_fit, mz_index)
mz_index$Value = mz_pred * 2

all_vote = extract_assigned_data(assigned_data, difference_cutoff = mz_index,
                                 difference_measure = "ObservedMZ", progress = TRUE)
