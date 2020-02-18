#' create lipids of interest
#'
#' Given a data.frame of lipid classifications, creates a data.frame which defines
#' the `isotopologue_EMF` as either of interest or not (1 / 0), which can
#' subsequently be used to define if the score for that EMF should be modified or not.
#'
#' @param lipid_classification data.frame of lipid classifications
#' @param lipid_weight what weighting to use for **lipids**
#' @param not_lipid_weight what weighting to use for **not_lipids** (metabolites)
#' @param other_weight what value to use for other things
#' @param all_emfs the full list of EMFs that we tried to classify
#'
#' @return data.frame
#' @seealso import_emf_classifications
#' @export
weight_lipid_classifications = function(lipid_classification, lipid_weight = 1,
                                     not_lipid_weight = 1, other_weight = 1,
                                     all_emfs = NULL){
  just_categories = unique(lipid_classification[, c("isotopologue_EMF", "Categories")])
  group_emfs = dplyr::group_by(just_categories, isotopologue_EMF)
  is_lipid = dplyr::summarise(group_emfs, lipid = sum(!(Categories %in% "not_lipid")) >= 1)
  is_lipid = dplyr::mutate(is_lipid, weight = dplyr::case_when(
    lipid ~ lipid_weight,
    !lipid ~ not_lipid_weight
  ))
  is_lipid$lipid = NULL

  if (!is.null(all_emfs)) {
    other_emfs = data.frame(isotopologue_EMF = all_emfs[!(all_emfs %in% is_lipid$isotopologue_EMF)],
                            weight = other_weight)

    out_lipids = rbind(is_lipid, other_emfs)
  } else {
    out_lipids = is_lipid
  }

}