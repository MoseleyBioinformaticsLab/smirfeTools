#' extracts IMF adducts
#'
#' Although technically all complete_EMFs are equivalent, they might not be
#' due to differences in the isotopologue adduct, for example Na vs K, 41K vs
#' 39K, H vs K vs Na. Therefore, we need to be able to save which adduct
#' belongs with each IMF.
#'
#' @param assignments a data.frame of assignments
#' @param imf which is the IMF column we want to keep
#' @param adduct which variable is holding the adduct information?
#'
#' @export
#' @return data.frame
keep_imf_adducts = function(assignments, imf = "complete_IMF", adduct = "adduct_IMF"){
  keep_data = assignments[assignments$Type %in% adduct, c(imf, "Assignment_Data")]
  names(keep_data)[2] = adduct
  keep_data
}