#' score and filter assignments
#'
#' Given assignments (either a filename or list), calculate scores, and filter
#' the assignments by M/Z
#'
#' @param assignments a character or list object
#' @param score_calculation how to calculate the score
#' @param score_weight multiplier to apply to the scores (default = 1)
#' @param filter_conditions how the data should be filtered
#' @param high_correlation what is considered high correlation with scan
#' @param emf_weight an optional data.frame of emf_weights
#' @param remove_highsd should the high standard deviation peaks be removed
#' @param remove_contaminants should peaks that were marked as contaminants be removed
#'
#' @seealso weight_lipid_classifications
#' @export
#' @return data.frame of assignments with scores added
score_filter_assignments = function(assignments, score_calculation = 1 - e_value,
                                    score_weight = 1, filter_conditions = ObservedMZ <= 1600,
                                    high_correlation = 0.5,
                                    emf_weight = NULL,
                                    remove_highsd = TRUE,
                                    remove_contaminants = TRUE){
  if (inherits(assignments, "character")) {
    assignments = readRDS(assignments)
  }

  sample_assignments = assignments$assignments
  sample_assignments = dplyr::mutate(sample_assignments, score = ({{score_calculation}}) * score_weight)

  sample_assignments = dplyr::filter(sample_assignments, {{filter_conditions}})

  sample_data = assignments$data
  # this tries to filter out assignments where the assignments are wholly due
  # to either peaks with high correlation with scan or peaks that have
  # high standard deviation in frequency space across scans. If all the peaks
  # in the assignment are like that, then the assignment is filtered out.
  # Any assignments where only some of the peaks are from those abberant conditions
  # are kept.
  if (all(c("ScanCorrelation", "HighScan", "HighSD") %in% names(sample_data)) && remove_highsd) {
    high_peaks = ((sample_data$ScanCorrelation >= high_correlation) & sample_data$HighScan) | sample_data$HighSD
    sample_data = sample_data[high_peaks, ]
    if (nrow(sample_data) > 0) {
      sample_assignments = sample_assignments[order(sample_assignments$complete_EMF), ]
      grouped_assignments = dplyr::group_by(sample_assignments, complete_EMF)
      summary_assignments = dplyr::summarise(grouped_assignments, n_peak = length(unique(PeakID)),
                                             n_high = sum(unique(PeakID) %in% sample_data$PeakID))
      reject_assignments = dplyr::filter(summary_assignments, n_peak == n_high)
      sample_assignments = dplyr::filter(sample_assignments, !(complete_EMF %in% reject_assignments$complete_EMF))
    }

  }
  if (("StandardContaminant" %in% names(sample_data)) && remove_contaminants) {
    standard_peaks = sample_data$PeakID[sample_data$StandardContaminant]
    sample_assignments = dplyr::mutate(sample_assignments, emf.size = paste0(complete_EMF, ".", clique_size))
    has_standard = sample_assignments$PeakID %in% standard_peaks
    filter_complete = unique(sample_assignments$complete_EMF[has_standard])
    grouped_possible_standard = dplyr::group_by(sample_assignments[sample_assignments$complete_EMF %in% filter_complete, ], emf.size)
    summary_hasstandard = dplyr::summarise(grouped_possible_standard, has_standard = any(unique(PeakID) %in% standard_peaks)) %>%
      dplyr::filter(has_standard)
    sample_assignments = dplyr::filter(sample_assignments, !(emf.size %in% summary_hasstandard$emf.size))
    sample_assignments$emf.size = NULL
  }

  # when a data.frame of *interesting* emfs is provided and a weighting, then the
  # scores for those emfs are further modified with a weighting factor
  if (!is.null(emf_weight)) {
    sample_assignments = dplyr::left_join(sample_assignments, emf_weight, by = "isotopologue_EMF")
    sample_assignments$weight[is.na(sample_assignments$weight)] = 1
    sample_assignments$score = sample_assignments$score * sample_assignments$weight
  }
  assignments$assignments = sample_assignments
  assignments
}

#' read SMIRFE assignments
#'
#' Given a SMIRFE JSON output file, read it in, along with other useful
#' information.
#'
#' @param smirfe_assignment the set of assignment results
#' @param assigned_only whether to return peaks with assignment only
#' @param .pb a progress bar object
#'
#' @importFrom jsonlite fromJSON
#'
#' @return list of tic, assignments, sample
#' @export
#'
read_smirfe_assignment <- function(smirfe_assignment, assigned_only = TRUE, .pb = NULL){
  log_memory()
  if (!is.null(.pb)) {
    knitrProgressBar::update_progress(.pb)
  }
  tmp_list <- jsonlite::fromJSON(smirfe_assignment, simplifyVector = FALSE)
  if (is.null(tmp_list$Sample)) {
    sample <- gsub(".json.output", "", basename(smirfe_assignment))
  } else {
    sample <- tmp_list$Sample
  }

  n_peaks = length(tmp_list$Peaks)

  if (!assigned_only) {
    to_extract = rep(TRUE, length(tmp_list$Peaks))
  } else {
    to_extract <- purrr::map_lgl(tmp_list$Peaks, function(x){length(x$Assignments) > 0})
  }
  n_extract = sum(to_extract)

  #tictoc::tic()
  extract_list = which(to_extract)
  extract_list = sample(extract_list, length(extract_list))
  peak_info <- internal_map$map_function(tmp_list$Peaks[extract_list], extract_peak_data)
  #tictoc::toc()
  peak_data <- purrr::map_dfr(peak_info, "peak_data")

  peak_assignments <- purrr::map_dfr(peak_info, "assignment")

  peak_data$Sample <- sample
  peak_data$Sample_Peak <- paste0(sample, "_", peak_data$PeakID)

  peak_assignments$Sample <- sample
  peak_assignments$Sample_Peak <- paste0(sample, "_", peak_assignments$PeakID)

  peak_assignments = dplyr::left_join(peak_assignments, peak_data[, c("PeakID", "ObservedMZ")], by = "PeakID")

  peak_id = paste0(sample, "_", unlist(tmp_list$ScanLevel$PeakID))
  keep_peaks = unique(peak_assignments$Sample_Peak)
  ideal_scanlevel = c("CorrectedLog10Height", "Log10Height",
                      "ObservedFrequency", "ObservedMZ")
  grab_scanlevel = ideal_scanlevel[ideal_scanlevel %in% names(tmp_list$ScanLevel)]
  scan_level_lists = purrr::map(tmp_list$ScanLevel[grab_scanlevel],
                                   function(in_matrix){
                                     tmp_list = purrr::map(in_matrix, ~ suppressWarnings(as.numeric(unlist(.x))))
                                     #tmp_matrix = do.call(rbind, tmp_list)
                                     names(tmp_list) = peak_id

                                     tmp_list[keep_peaks]
                                   })
  log_memory()

  list(tic = tmp_list$TIC,
       assignments = peak_assignments,
       data = peak_data,
       scan_level = scan_level_lists,
       sample = sample,
       peaks = data.frame(sample = sample,
                          n_json = n_peaks,
                          n_extracted = n_extract,
                          stringsAsFactors = FALSE))
}

get_tic <- function(assigned_data){
  tic <- purrr::map(assigned_data, "tic")
  tic = unlist(purrr::map(tic, as.double))
  names(tic) <- purrr::map_chr(assigned_data, "sample")
  tic
}

#' Remove Only Labeled EMFs
#'
#' Removes EMFs that *only* contain labeled IMFs.
#'
#' When working with metabolomics data where labeled precursors are used
#' and there is a reasonable expectation that label is incorporated,
#' it is highly unlikely that there are cases where the un-labeled
#' IMF among a set of EMFs will not be observed.
#' However, this is not baked into the
#' SMIRFE assignment algorithm, so we provide this functionality here.
#'
#' Note that this should be run **after** `read_smirfe_assignment` and
#' **before** `extract_assigned_data`.
#'
#' @param assignment_data a data.frame of assignments for peaks
#' @param remove_s logical, should EMFs with Sulfur also be removed?
#'
#' @export
#' @return data.frame
remove_only_labeled_emfs = function(assignment_data, remove_s = TRUE){
  if (remove_s) {
    s_adducts = grepl("S", assignment_data$complete_EMF)
    assignment_data = assignment_data[!s_adducts, ]
  }

  lbl_counts = assignment_data[assignment_data$Type %in% "lbl.count", ]
  split_emf = split(lbl_counts, lbl_counts$complete_EMF)

  keep_emf = purrr::map_lgl(split_emf, function(in_adduct){
    if (sum(in_adduct$Assignment_Data == 0) == 0){
      return(FALSE)
    } else {
      return(TRUE)
    }
  })

  valid_emf = names(keep_emf)[keep_emf]

  assignment_data = assignment_data[assignment_data$complete_EMF %in% valid_emf, ]
  assignment_data
}


#' extract data
#'
#' Extract the data from a complete set of assignments.
#'
#' @param assigned_data a list of assignments created from `read_smirfe_assignments`
#' @param emf_classifications data.frame of isotopologue_EMF classifications
#' @param sample_peak which variable holds the sample peak
#' @param imf which variable holds the IMF information
#' @param e_value the e-values
#' @param evalue_cutoff only consider e_values below this cutoff
#' @param emf which variable holds the EMF information
#' @param sample which variable holds the sample
#' @param observed_mz which variable holds the observed M/Z
#' @param observe_frequency which variable holds the observed frequency data
#' @param assigned_mz which variable holds the assigned M/Z
#' @param height which variable holds the height / intensity
#' @param other_cols which other columns should be kept?
#' @param chosen_keep_ratio what is the ratio of max to keep chosen EMFs?
#' @param difference_cutoff what is the maximum difference that peaks in an IMF should have?
#' @param difference_measure which set of values from the data should be used to check differences?
#' @param use_scan_level should scan level data be used for estimating SDs?
#' @param progress should the progress be shown?
#'
#' @export
#'
#' @importFrom purrr map map_df map_int
#' @importFrom dplyr left_join
#'
#' @return list of matrices and a data.frame
extract_assigned_data <- function(assigned_data,
                                  emf_classifications = NULL,
                                  sample_peak = "Sample_Peak",
                                  imf = "complete_IMF",
                                  e_value = "e_value",
                                  evalue_cutoff = 0.5,
                                  emf = "complete_EMF",
                                  sample = "Sample",
                                  mz_diff = "mass_error",
                                  data_col = "Assignment_Data",
                                  numeric_values = c("e_value", "mass_error", "NAP", "lbl.count", "clique_size"),
                                  observed_mz = "ObservedMZ",
                                  observed_frequency = "ObservedFrequency",
                                  height = "Height",
                                  remove_elements = "S",
                                  chosen_keep_ratio = 0.9,
                                  difference_cutoff = NULL,
                                  difference_measure = "ObservedFrequency",
                                  use_scan_level = TRUE,
                                  progress = TRUE){

  if (is.null(assigned_data[[1]]$assignments$score)) {
    stop("No 'scores' are available, do you need to run `calculate_assignment_scores`?")
  }
  if (is.null(difference_cutoff)) {
    warning("difference_cutoff is NULL, no checking will be done on IMFs!")

    if ("ObservedFrequency" %in% difference_measure) {
      warning("difference_cutoff set to 1 frequency point difference!")
    } else {
      warning("difference_cutoff will be set to 1ppm of M/Z!")
    }
    difference_cutoff = NA
    names(difference_cutoff) = difference_measure
  }

  if (!use_scan_level) {
    warning("use_scan_level is FALSE, not using scan information!")
    scan_level_location = NULL
  } else {
    scan_level_location = purrr::map(assigned_data, ~ .x$scan_level[[difference_measure]])
    scan_level_names = unlist(purrr::map(scan_level_location, ~ names(.x)))
    scan_level_location = unlist(scan_level_location, recursive = FALSE, use.names = FALSE)
    names(scan_level_location) = scan_level_names
  }

  start_time = Sys.time()
  if (progress) {
    message("Generating EMF cliques from each sample ...")
  }
  log_message("Generating EMF cliques from each sample ...")

  names(assigned_data) = NULL

  within_sample_emfs = internal_map$map_function(assigned_data, function(.x){
    progress_msg = paste0("Sample: ", .x$sample)
    if (progress) {
      message(progress_msg)
    }
    log_message(progress_msg)
    tmp_assign = dplyr::filter(.x$assignments, !grepl(remove_elements, complete_EMF))
    get_sample_emfs(tmp_assign, .x$sample, evalue_cutoff = evalue_cutoff, emf_classifications = emf_classifications)
  })


  all_gemf_emf_mapping = internal_map$map_function(within_sample_emfs, function(x){
    purrr::map_df(x, ~ unique(dplyr::select(.x, grouped_EMF, complete_EMF)))
  })

  all_gemf_emf_mapping = do.call(rbind, all_gemf_emf_mapping)

  if (progress) {
    message("Creating pseudo EMFs across cliques ...")
  }
  log_message("Creating pseudo EMFs across cliques ...")
  log_memory()
  sudo_emf_list = create_sudo_emfs(all_gemf_emf_mapping)
  # next things:
  all_gemfs = unlist(within_sample_emfs, recursive = FALSE, use.names = FALSE)
  names(all_gemfs) = purrr::map_chr(all_gemfs, ~ .x$grouped_EMF[1])

  n_gemf = length(all_gemfs)

  if (progress) {
    message("Choosing EMFs by voting ...")
  }
  log_message("Choosing EMFs by voting ...")

  get_measure = function(data, measure){
    data %>% dplyr::select(.data[[measure]], Sample_Peak)
  }
  peak_location = purrr::map_df(assigned_data, function(in_data){
    get_measure(in_data$data, difference_measure)
  })
  match_name = names(peak_location) %in% difference_measure
  names(peak_location)[match_name] = "Value"

  sudo_emf_list = sudo_emf_list[sample(length(sudo_emf_list), length(sudo_emf_list))]

  n_sudoemf = length(sudo_emf_list)

  chosen_emfs = internal_map$map_function(sudo_emf_list, function(.x){
    #message(.x$sudo_EMF[1])
    choose_emf(all_gemfs[unique(.x$grouped_EMF)], scan_level_location, peak_location, difference_cutoff, chosen_keep_ratio, .x$sudo_EMF[1])
  })
  names(chosen_emfs) = names(sudo_emf_list)
  null_chosen = purrr::map_lgl(chosen_emfs, ~ nrow(.x) == 0)
  chosen_emfs = chosen_emfs[!null_chosen]

  n_chosen = length(chosen_emfs)

  # debugging version
  # chosen_emfs = purrr::map(seq_along(sudo_emf_list), function(.x){
  #   message(.x)
  #   choose_emf(all_gemfs[unique(sudo_emf_list[[.x]]$grouped_EMF)], scan_level_location, peak_location, difference_cutoff, chosen_keep_ratio)
  # })


  if (progress) {
    message("Merging chosen EMFs ...")
  }
  log_message("Merging chosen EMFs ...")
  merged_chosen_emfs = merge_duplicate_semfs(chosen_emfs, all_gemfs, scan_level_location, peak_location, difference_cutoff, chosen_keep_ratio)
  # next is to actually extract the right data. But up to here, everything appears OK.
  #
  null_chosen2 = purrr::map_lgl(merged_chosen_emfs, ~ nrow(.x) == 0)
  merged_chosen_emfs = merged_chosen_emfs[!null_chosen2]

  n_merged = length(merged_chosen_emfs)
  log_memory()
  if (progress) {
    message("Extracting EMF matrices ...")
  }
  log_message("Extracting EMF matrices ...")
  extracted_emfs = extract_emfs(merged_chosen_emfs)

  all_peaks = purrr::map_df(merged_chosen_emfs, ~ unique(dplyr::select(.x, Sample, Sample_Peak)))
  peak_data = purrr::map_df(assigned_data, ~ .x$data)
  peak_data = dplyr::filter(peak_data, Sample_Peak %in% all_peaks$Sample_Peak)

  extracted_location_intensity = add_location_intensity(extracted_emfs, peak_data,
                                                        ObservedMZ, Height)

  all_assignments = purrr::map_df(merged_chosen_emfs, ~ .x)

  stop_time = Sys.time()
  diff_time = difftime(stop_time, start_time)
  time_message = paste0("Done in ", diff_time, " ", attr(diff_time, "units"))
  if (progress) {
    message(time_message)
  }
  log_message(time_message)
  n_emfs = data.frame(type = c(
    "grouped",
    "sudo",
    "chosen",
    "merged"
  ),
  number = c(
    n_gemf,
    n_sudoemf,
    n_chosen,
    n_merged
  ),
  stringsAsFactors = FALSE)
  return(list(emfs = extracted_location_intensity,
              emf_info = all_assignments,
              tic = get_tic(assigned_data),
              n_emfs = n_emfs)
  )
}

filter_peak_imfs <- function(peak_assignments, imf = "complete_IMF", e_value = "e_value"){
  e_value <- rlang::enquo(e_value)
  #imf <- rlang::enquo(imf)

  if (length(unique(peak_assignments[[imf]])) == 1) {
    return(peak_assignments)
  } else {
    just_evalues <- dplyr::filter(peak_assignments, Type %in% !!e_value)
    just_evalues$Assignment_Data <- as.numeric(just_evalues$Assignment_Data)
    just_evalues <- just_evalues[order(just_evalues$Assignment_Data), ]
    log_evalues <- -1*log10(just_evalues$Assignment_Data)
    evalue_diffs <- log_evalues - dplyr::lead(log_evalues)
    evalue_diffs[is.infinite(evalue_diffs)] <- 1e6
    evalue_diffs[is.na(evalue_diffs)] <- 0

    if (sum(evalue_diffs >= 3) > 0) {
      keep_imf <- just_evalues[[imf]][seq(1, which(evalue_diffs >= 3)[1])]
    } else {
      keep_imf <- just_evalues[[imf]]
    }
    out_peaks <- peak_assignments[peak_assignments[[imf]]%in% keep_imf, ]
    return(out_peaks)
  }
}

extract_single <- function(assignment_data, peaks, use_var = NULL, key = NULL, value = NULL, sample_peak = "Sample_Peak", sample = "Sample"){

  if (is.null(use_var)) {
    stop("The variable to extract (use_var) must be defined!")
  }

  tmp_data <- assignment_data[(assignment_data[[key]] %in% use_var) & (assignment_data[[sample_peak]] %in% peaks), ]
  if (!is.na(as.numeric(tmp_data[[value]][1]))) {
    tmp_data[[value]] <- as.numeric(tmp_data[[value]])
  }

  out_data <- tmp_data[[value]]
  names(out_data) <- tmp_data[[sample]]

  out_data

}

extract_multiple <- function(assignment_data, use_var = NULL, sample = "sample"){

  if (is.null(use_var)) {
    stop("The variable to extract (use_var) must be defined!")
  }

  tmp_data <- split(assignment_data[, use_var], assignment_data[, sample])
  tmp_data <- purrr::map(tmp_data, unique)
  tmp_data
}

#' choose a peak
#'
#' When there are multiple peaks in a sample, we need to choose one based on
#' the chosen representative IMFs. See **Details** for more information on how
#' that decision is made.
#'
#' @param mz_diff_data the data.frame with mz_errors
#' @param imf which variable are the IMFs listed
#' @param sample_peak which variable contains the sample peak information
#' @param sample which variable has the sample information
#' @param data_col which variable has the actual mz mass error
#'
#' @details If there are multiple peaks for each sample, picks the one with
#' the smallest M/Z error with the provided isotopologues.
#'
#' @export
#'
#' @return data.frame Peaks and Samples
choose_single_peak <- function(mz_diff_data, imf = "IMF",
                               sample_peak = "Sample_Peak", sample = "Sample",
                               data_col = "Assignment_Data"){
  split_sample <- split(mz_diff_data, mz_diff_data[, sample])

  keep_peaks <- purrr::map_chr(split_sample, function(in_sample){
    peak_index <- which.min(in_sample[[data_col]])
    in_sample[[sample_peak]][peak_index]
  })
  keep_peaks
}


#' decide on IMFs
#'
#' given a set of assignment data, decide which set of IMFs should
#' be reported.
#'
#' @param assignment_data data.frame of assignment data
#' @param sample_peak what defines a sample_peak combination
#' @param imf which is the imf variable
#' @param sample what is the sample
#' @param e_value which is the e-value column?
#' @param min_e_value what is the minimum e-value to use?
#'
choose_imfs <- function(assignment_data, sample_peak = "Sample_Peak", imf = "complete_IMF",
sample = "Sample", e_value = "e_value", data_col = "Assignment_Data"){

  # easy cases, there is only 1 IMF
  if (length(unique(assignment_data[[imf]])) == 1) {
    return(assignment_data[, imf])
  }

  assignment_e_values <- assignment_data[assignment_data$Type %in% e_value, ]
  assignment_e_values[[data_col]] <- as.numeric(assignment_e_values[[data_col]])
  # then there are multiple IMFs to choose from
  assignment_single_e_value <- pick_single_evalue(assignment_e_values,
                                                  sample_peak = sample_peak,
                                                  imf = imf,
                                                  sample = sample,
                                                  data_col = data_col)

  split_assignment_data <- split(assignment_single_e_value, assignment_single_e_value[, imf])
  n_sample <- length(unique(assignment_single_e_value[, sample]))

  imf_stats <- purrr::map_df(split_assignment_data, function(in_imf){
    data.frame(imf = unique(in_imf[, imf]),
               ratio = nrow(in_imf) / n_sample,
               evalue = median_e_values(in_imf[, data_col]),
               stringsAsFactors = FALSE)
  })

  # logic:
  # case 1: everything is possible:
  #   * all ratios of max ratio to everything else is < 0.2
  #   * all differences of min(log10(evalue)) to everything else is < 3, so
  #   everything is valid.
  # case 2: multiple valid IMFs:
  #   * max ratio to next ratio is > 0.2
  #   * there is a difference in log10(evalue)'s > 3


  max_ratio_value <- max(imf_stats$ratio)
  ratio_diffs <- abs(max_ratio_value - imf_stats$ratio)

  out_ratio <- !(ratio_diffs > 0.2)


  if (all(imf_stats$evalue == 0)) {
    out_evalue <- rep(TRUE, nrow(imf_stats))
  } else {
    min_evalue <- min(imf_stats$evalue)
    if (min_evalue == 0) {
      min_evalue <- 1
    }
    imf_stats[imf_stats$evalue == 0, "evalue"] <- 1

    evalue_diffs <- abs(-1*log10(min_evalue) - -1*log10(imf_stats$evalue))
    out_evalue <- !(evalue_diffs > 3)
  }


  choose <- out_ratio & out_evalue

  # if we've got nothing, try the next highest ratio
  if (sum(choose) == 0) {
    max_ratio_value <- max(imf_stats[imf_stats$ratio != max_ratio_value, "ratio"])
    ratio_diffs <- abs(max_ratio_value - imf_stats$ratio)

    out_ratio <- !(ratio_diffs > 0.2)
    choose <- out_ratio & out_evalue
  }

  # if still nothing, do the | (technically, this should always get something)
  if (sum(choose) == 0) {
    choose <- out_ratio | out_evalue
  }

  # but just in case, return everything if that is still nothing
  if (sum(choose) == 0) {
    choose <- rep(TRUE, nrow(imf_stats))
  }

  return(imf_stats[choose, "imf"])

}

pick_single_evalue <- function(assignment_e_values, sample_peak = "Sample_Peak", imf = "complete_IMF",
                               sample = "Sample", data_col = "Assignment_Data"){
  split_imf <- split(assignment_e_values, paste0(assignment_e_values[, imf], assignment_e_values[, sample]))

  assignment_evalues_single <- purrr::map_df(split_imf, function(x){
    nrow_unique <- nrow(unique(x[, c(imf, data_col)]))
    if (nrow_unique > 1) {
      return(x[which.min(x[, data_col]), ])
    } else {
      return(x[1, ])
    }
  })
}


median_e_values <- function(e_values){
  # handle na
  n_value <- length(e_values)
  if (sum(is.na(e_values)) == n_value) {
    return(1)
  } else {
    e_values <- e_values[!is.na(e_values)]
  }
  #e_values[e_values == 0] <- min_e_value
  if (n_value < 15) {
    min_value <- min(e_values)
    return(min_value)
  } else {
    return(median(e_values))
  }
}


#' create pseudo peak
#'
#' Given a master table of assignments, creates a set of "pseudo" peaks
#' that each contain the set of `IMF`s and the `sample_peak` identifiers
#' to pull out of a master table of peak assignments.
#'
#' @param in_assignments the table of all assignments from all samples
#' @param sample_peak which variable holds the sample peak
#' @param imf which variable holds the IMF information
#'
#' @importFrom purrr map_lgl
#'
#' @export
#' @return list of pseudo peaks
create_sudo_peaks <- function(in_assignments, sample_peak = "Sample_Peak", imf = "complete_IMF"){
  in_assignments$grabbed <- FALSE

  peak_index <- 1

  # makes it possible to reuse the variable names in other places where we
  # need matching between this list and the assignments data.frame
  mock_list <- vector("list", length = 2)
  names(mock_list) <- c(sample_peak, imf)

  # maximum amount of peaks, likely much less, but this sets an upper bound for
  # us.
  sudo_peaks <- vector("list", length(unique(in_assignments[!in_assignments$grabbed, sample_peak])))

  # while there is something grab, go in and get the first peak. We then grab
  # the IMFs for that sample peak, then all the sample peaks for those IMFs,
  # continuing until the set of IMFs and sample peaks no longer change.
  #
  # This should be equivalent to doing set intersections + union across a list
  # of IMFs for sample peaks, but this is much, much faster than checking all
  # of the sample peaks
  while (sum(!in_assignments$grabbed) > 0) {
    tmp_peak <- unique(in_assignments[!in_assignments$grabbed, sample_peak])[1]
    tmp_imf <- in_assignments[in_assignments[, sample_peak] %in% tmp_peak, imf][1]

    all_peak_imf <- unique(in_assignments[in_assignments[, imf] %in% tmp_imf, sample_peak])
    all_imf_peak <- unique(in_assignments[in_assignments[, sample_peak] %in% all_peak_imf, imf])

    match_peaks <- FALSE

    n_iter <- 1
    while (!match_peaks) {
      all_peak_imf_2 <- unique(in_assignments[in_assignments[, imf] %in% all_imf_peak, sample_peak])
      all_imf_peak_2 <- unique(in_assignments[in_assignments[, sample_peak] %in% all_peak_imf_2, imf])

      match_peaks <- isTRUE(all.equal(all_peak_imf_2, all_peak_imf)) &&
        isTRUE(all.equal(all_imf_peak_2, all_imf_peak))

      all_peak_imf <- all_peak_imf_2
      all_imf_peak <- all_imf_peak_2
      n_iter <- n_iter + 1
    }
    #print(n_iter)
    tmp_list <- mock_list
    tmp_list[[sample_peak]] <- unique(all_peak_imf)
    tmp_list[[imf]] <- unique(all_imf_peak)
    sudo_peaks[[peak_index]] <- tmp_list

    in_assignments[in_assignments[, imf] %in% all_imf_peak, "grabbed"] <- TRUE
    peak_index <- peak_index + 1

  }
  null_sudos <- purrr::map_lgl(sudo_peaks, is.null)
  sudo_peaks <- sudo_peaks[!null_sudos]
  sudo_peaks
}


#' peak mz, height, assignments
#'
#' Extracts desired peak variables and assignment information from a peak list,
#' the peak list is typically read in from a JSON output file of assignments
#' generated by SMIRFE.
#'
#' @param peak_list the peak list from reading some JSON
#' @param variables which variables to extract (default: ObservedMZ, Height)
#' @param summaries which summary variables to extract (default: Mean)
#' @param include_non_primary include non-primary assignments (default: TRUE)
#' @param keep_peaks keep "both" original and assigned peaks, "assignments" only, or original "peaks"
#'
#' @export
#' @return data.frame
#'
get_assigned_peak_info <- function(peak_list, variables = c("ObservedMZ", "Height"),
                                   summaries = "Mean", include_non_primary = TRUE,
                                   keep_peaks = "assignments"){

  peak_df_1 <- peak_list_2_df(peak_list, summary_values = variables, summary_types = summaries)
  peak_df_2 <- extract_peak_characteristics(peak_df_1, extract_variables = variables, extract_summary = summaries)

  peak_assignments <- peak_assignments_2_df(peak_list, include_non_primary = include_non_primary)

  merge_peak_characteristics_assignments(peak_df_2, peak_assignments, keep_only = keep_peaks)
}

peak_list_2_df <- function(peak_list,
                           summary_values = c("ObservedMZ", "Height", "Area", "NormalizedArea"),
                           summary_types = c("Mean", "Median", "SD", "RSD", "ModelSD")) {
  n_peaks <- length(peak_list)
  peak_index <- seq(1, n_peaks)

  purrr::map_df(peak_index, function(ipeak){
    #print(ipeak)
    single_peak <- peak_list[[ipeak]]

    if (is.null(single_peak$PeakID)) {
      peak_id <- ipeak
    } else {
      peak_id <- single_peak$PeakID
    }
    # if (is.na(single_peak$Sample) & !is.null(single_peak$Scan)) {
    #   scan <- as.character(single_peak$Scan)
    # } else {
    #   scan <- single_peak$Sample
    # }
    nscan <- single_peak$NScan

    tmp_df <- purrr::map_df(summary_values, function(in_value){
      use_summary <- intersect(names(single_peak[[in_value]]), summary_types)
      summary_df <- purrr::map_df(use_summary, function(in_summary){
        data.frame(value = single_peak[[in_value]][[in_summary]],
                   summary = in_summary, stringsAsFactors = FALSE)
      })
      summary_df$extracted <- in_value
      summary_df
    })
    tmp_df$peak_id <- peak_id
    tmp_df$scan <- scan
    tmp_df$nscan <- nscan

    tmp_df
  })
}

extract_assignments <- function(assignment_list){
  if (!is.null(assignment_list$lbl)) {
    if (length(assignment_list$lbl) == 1) {
      assignment_list$lbl = NULL
    } else if (length(assignment_list$lbl) == 2) {
      names(assignment_list$lbl) <- c("type", "count")
    }

  }

  not_null_entries = !purrr::map_lgl(assignment_list, is.null)
  assignment_list = assignment_list[not_null_entries]
  assignment_df <- as.data.frame(assignment_list, stringsAsFactors = FALSE)

}

extract_peak_data <- function(single_peak_list){
  "!DEBUG PeakID = `single_peak_list$PeakID`"
  assignments <- single_peak_list$Assignments
  single_peak_list$Assignments <- NULL
  if (!is.null(single_peak_list$type)) {
    single_peak_list$type <- NULL
  }

  single_names <- names(single_peak_list)
  keep_names <- single_names[!(single_names %in% "PeakID")]


  single_df <- as.data.frame(single_peak_list, stringsAsFactors = FALSE)

  if (length(assignments) != 0) {
    peak_assignment <- suppressWarnings(purrr::map_df(assignments, extract_assignments))
    peak_assignment$PeakID <- single_df$PeakID[1]
  } else {
    peak_assignment <- list()
  }

  return(list(peak_data = single_df, assignment = peak_assignment))

}

extract_peak_characteristics <- function(peak_df, extract_variables = c("ObservedMZ", "Height"), extract_summary = c("Mean")){
  expanded_combinations <- expand.grid(extract_variables, extract_summary)
  expected_variables <- paste0(expanded_combinations[,1], ".", expanded_combinations[,2])

  peak_df <- peak_df[(peak_df$extracted %in% extract_variables) & (peak_df$summary %in% extract_summary), ]
  peak_df$extract_variable <- paste0(peak_df$extracted, ".", peak_df$summary)

  split_by_id <- split(peak_df, peak_df$peak_id)
  purrr::map_df(split_by_id, function(x){
    tmp_vector <- x$value
    names(tmp_vector) <- x$extract_variable
    out_vector <- tmp_vector[expected_variables]
    names(out_vector) <- expected_variables
    out_frame <- as.data.frame(as.list(out_vector), stringsAsFactors = FALSE)
    out_frame$peak_id <- x$peak_id[1]
    out_frame$nscan <- x$nscan[1]
    out_frame$scan <- x$scan[1]
    out_frame
  })
}

merge_peak_characteristics_assignments <- function(extracted_characteristics, assignment_df, keep_only = "assignments"){
  join_type <- switch(keep_only,
                      assignments = dplyr::right_join,
                      peaks = dplyr::left_join,
                      both = dplyr::full_join)

  out_df <- join_type(extracted_characteristics, assignment_df, by = c("peak_id"))

  out_df
}


peak_assignments_2_df <- function(peak_list, include_non_primary = TRUE) {
  n_peaks <- length(peak_list)
  peak_index <- seq(1, n_peaks)

  assignments <- purrr::map(peak_index, function(ipeak){
    #print(ipeak)
    single_peak <- peak_list[[ipeak]]

    if (is.null(single_peak$Peak)) {
      peak_id <- ipeak
    } else {
      peak_id <- single_peak$Peak
    }

    if (length(single_peak$Assignments) != 0) {
      assign_df <- assignment_list_2_df(single_peak$Assignments)
      assign_df$peak_id <- peak_id
    } else {
      assign_df <- NULL
    }
    assign_df
  })

  null_assignments <- purrr::map_lgl(assignments, is.null)
  assignments <- assignments[!null_assignments]
  assignments <- purrr::map_df(assignments, function(x){x})
  assignments <- cbind(assignments, split_adduct(assignments[["EMF"]]))
  assignments
}


assignment_list_2_df <- function(in_assignment, include_non_primary = TRUE){
  if (!is.null(in_assignment[["Type"]])) {
    assign_types <- purrr::map_chr(in_assignment, "Type")

    if (!include_non_primary) {
      in_assignment <- in_assignment[assign_types %in% "Primary"]
    }

    in_assignment <- purrr::map(in_assignment, replace_null)

  }

  assign_df <- purrr::map_df(in_assignment, as.data.frame, stringsAsFactors = FALSE)

}

split_adduct <- function(emf, split_chr = "_"){
  emf_split <- strsplit(emf, split = split_chr, fixed = TRUE)
  purrr::map_df(emf_split, function(x){
    data.frame(Adduct = x[1], EMF2 = x[2], stringsAsFactors = FALSE)
  })
}

replace_null <- function(in_list, replace_value = NA) {
  purrr::map(in_list, function(x){
    if (is.null(x)) {
      replace_value
    } else {
      x
    }
  })
}


#' Complete EMF mapping
#'
#' Everything is reported by `complete_EMF` or `complete_IMF`, which includes
#' the adducts, to make sure we have unique peaks. However, for some other
#' work, we actually want the **non-complete** EMF and IMF. This function
#' generates the mapping of `complete_EMF` to `isotopologue_EMF / IMF` for
#' other uses.
#'
#' @param assigned_data lists of assigned data
#' @param remove_S should formulas containing "S" be removed
#' @param progress show a progress bar for the main worker
#'
#' @return data.frame
#' @export
complete_emf_mappings = function(assigned_data, remove_S = TRUE, progress = TRUE){
  all_assignments = internal_map$map_function(assigned_data, function(x){
    tmp_assign = x$assignments
    tmp_assign = dplyr::filter(tmp_assign,
                  Type %in% c("isotopologue_EMF",
                              "isotopologue_IMF"))
    if (remove_S) {
      tmp_assign = dplyr::filter(tmp_assign,
                                 !(grepl("S", complete_EMF)))
    }
    tmp_assign = dplyr::select(tmp_assign, -PeakID,
                  -Sample, -Sample_Peak)
    tmp_assign = unique(tmp_assign)
    tmp_assign
  })
  all_assignments = do.call(rbind, all_assignments)
  rownames(all_assignments) = NULL

  all_assignments = unique(all_assignments)
  all_assignments$emf_imf = paste0(all_assignments$complete_EMF, ".", all_assignments$complete_IMF)
  split_assignments = split(all_assignments, all_assignments$emf_imf)
  #split_assignments = unique(split_assignments)

  if (progress) {
    pb = knitrProgressBar::progress_estimated(length(split_assignments))
  } else {
    pb = NULL
  }

  spread_assignments = purrr::map_dfr(split_assignments, function(in_assign, .pb){
    knitrProgressBar::update_progress(.pb)
    if (nrow(in_assign) == 2) {
      return(tidyr::spread(in_assign, Type, Assignment_Data))
    } else {
      assign_index = seq(1, nrow(in_assign)-1, 2)
      out_assign = purrr::map_df(assign_index, function(in_index){
        tidyr::spread(in_assign[in_index:(in_index + 1), ],
                      Type, Assignment_Data)

      })
      return(out_assign)
    }
  }, pb)
  spread_assignments
}
