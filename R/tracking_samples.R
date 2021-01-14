#' remove samples
#'
#' Provides ways for the user to note in a sample data.frame which samples were
#' kept or removed, and for which reasons.
#'
#' @param sample_df the sample data.frame to modify
#' @param logical_condition which samples are to be removed
#' @param reason why should the samples be removed
#'
#' @export
#' @return data.frame
remove_samples = function(sample_df,
                          logical_condition = NULL,
                          reason = NULL){

  if (is.null(logical_condition)) {
    warning("No filtering applied!")
  }

  if (is.null(reason)) {
    stop("The 'reason' cannot be NULL, please supply a character string, eg 'Something was wrong with these'")
  }

  if (nchar(reason) == 0) {
    stop("The 'reason' can't be an empty string (''), please supply some information about why the samples are filtered out")
  }

  if (length(logical_condition) != nrow(sample_df)) {
    stop("The 'logical_condition' has to be the same length as the 'sample_df' supplied.")
  }

  if (!("removed" %in% names(sample_df))) {
    sample_df$removed = "kept"
  }

  sample_df$removed[which(logical_condition)] = reason
  sample_df
}

#' create sample tracking
#'
#' Create a data.frame for sample tracking and filtering.
#'
#' @param sample_list list of samples
#' @param contains_paths do the samples contain full paths
#'
#' @export
#' @return data.frame
create_sample_tracking = function(sample_list){
  sample_df = data.frame(sample_id = basename(sample_list),
                         sample_path = sample_list)
  sample_df
}