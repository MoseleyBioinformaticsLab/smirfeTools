#' calculate mean fold changes
#'
#' Provided values of intensities and groups of samples, calculates log2 fold-changes
#' between the mean values across the groups of samples.
#'
#' @param intensities matrix of intensities, rows are features, columns are samples
#' @param grouped_samples list of samples
#' @param add_min a minimum value to add to 0 so that log doesn't complain
#'
#' @importFrom purrr map_dfc
#'
#' @return data.frame
#' @export
calculate_mean_fc = function(intensities, grouped_samples, add_min = 1e-12){
  added_int = intensities + add_min

  group_means = purrr::map_dfc(grouped_samples, ~ log2(rowMeans(added_int[, .x])))
  names(group_means) = names(grouped_samples)

  lfc = data.frame(lfc = group_means[, 2] - group_means[, 1])
  names(lfc) = paste0("Log2.FC.", names(grouped_samples)[2], ".", names(grouped_samples)[1])

  out_fc = cbind(group_means, lfc)
  out_fc$feature = rownames(intensities)
  out_fc
}