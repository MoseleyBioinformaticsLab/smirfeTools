#' read SMIRFE assignments
#'
#' Given a SMIRFE JSON output file, read it in, along with other useful
#' information.
#'
#' @param smirfe_assignment the set of assignment results
#' @param .pb a progress bar object
#'
#' @importFrom jsonlite fromJSON
#'
#' @return list of tic, assignments, sample
#' @export
#'
read_smirfe_assignment <- function(smirfe_assignment, .pb = NULL){
  if (!is.null(.pb)) {
    knitrProgressBar::update_progress(.pb)
  }
  tmp_list <- jsonlite::fromJSON(smirfe_assignment, simplifyVector = FALSE)
  sample <- gsub(".json.output", "", basename(smirfe_assignment))

  assignments <- get_assigned_peak_info(tmp_list$Peaks)
  assignments$sample <- sample
  assignments$sample_peak <- paste0(assignments$sample, "_", assignments$peak_id)
  list(tic = tmp_list$TotalIntensity$Value,
       assignments = assignments,
       sample = sample)
}

get_tic <- function(assigned_data){
  tic <- purrr::map_dbl(assigned_data, "tic")
  names(tic) <- purrr::map_chr(assigned_data, "sample")
  tic
}

#' extract data
#'
#' Extract the data from a complete set of assignments.
#'
#' @param assigned_data a list of assignments created from `read_smirfe_assignments`
#' @param remove_secondary remove the seconardy assignments? (default is TRUE)
#' @param sample_peak which variable holds the sample peak
#' @param imf which variable holds the IMF information
#' @param e_value the e-values
#' @param emf which variable holds the EMF information
#' @param sample which variable holds the sample
#' @param observed_mz which variable holds the observed M/Z
#' @param assigned_mz which variable holds the assigned M/Z
#' @param height which variable holds the height / intensity
#' @param other_cols which other columns should be kept?
#' @param progress should the progress be shown?
#'
#' @export
#'
#' @importFrom purrr map map_df map_int
#' @importFrom dplyr left_join
#'
#' @return list of matrices and a data.frame
extract_assigned_data <- function(assigned_data, remove_seconary = TRUE,
                                  sample_peak = "sample_peak",
                                  imf = "IMF", e_value = "IMF.E.Value",
                                  emf = "EMF", sample = "sample",
                                  observed_mz = "ObservedMZ.Mean",
                                  assigned_mz = "Assigned.M.Z",
                                  height = "Height.Mean",
                                  other_cols = c("Adduct_IMF", "Adduct", "EMF2"),
                                  progress = TRUE){
  all_assignments <- purrr::map_df(assigned_data, "assignments")

  if (remove_seconary) {
    all_assignments <- all_assignments[all_assignments$Type %in% "Primary", ]
  }


  if (progress) message("Creating pseudo peaks ...")

  sudo_start <- Sys.time()
  sudo_peaks <- create_sudo_peaks(all_assignments, sample_peak = sample_peak, imf = imf)
  sudo_end <- Sys.time()

  if (progress) difftime(sudo_end, sudo_start)

  all_e_value <- all_assignments[, e_value]
  min_e_value <- min(all_e_value[all_e_value > 0], na.rm = TRUE)
  min_e_value <- min_e_value * 0.1

  all_samples <- unique(all_assignments[, sample])
  peak_matrix <- matrix(NA, nrow = length(sudo_peaks), ncol = length(all_samples))
  colnames(peak_matrix) <- all_samples
  rownames(peak_matrix) <- paste0("X", seq(1, length(sudo_peaks)))

  height_matrix <- peak_matrix
  mz_matrix <- peak_matrix
  imf_matrix <- matrix(vector("list", 1), nrow = length(sudo_peaks), ncol = length(all_samples))
  colnames(imf_matrix) <- all_samples
  rownames(imf_matrix) <- rownames(peak_matrix)
  emf_matrix <- imf_matrix

  peak_2_imf <- vector("list", length(sudo_peaks))
  names(peak_2_imf) <- rownames(peak_matrix)
  peak_2_peak <- peak_2_imf

  if (progress) message("Extracting data ...")

  if (progress && require(knitrProgressBar)) {
    pb <- knitrProgressBar::progress_estimated(length(sudo_peaks))
    do_update <- TRUE
  } else {
    do_update <- FALSE
  }

  for (ipeak in seq_along(sudo_peaks)) {
    if (do_update) {
      knitrProgressBar::update_progress(pb)
    }


    tmp_assignment <- all_assignments[all_assignments[, imf] %in% sudo_peaks[[ipeak]][[imf]], ]
    use_imfs <- choose_imfs(tmp_assignment, sample_peak = sample_peak,
                            imf = imf, sample = sample, e_value = e_value,
                            min_e_value = min_e_value)
    peak_2_imf[[ipeak]] <- use_imfs

    imf_mz <- unique(tmp_assignment[tmp_assignment[, imf] %in% use_imfs, assigned_mz])
    peak_data <- choose_single_peak(tmp_assignment, use_imfs, imf_mz, imf = imf,
                                    sample_peak = sample_peak, sample = sample)

    use_peaks <- unique(peak_data[, sample_peak])
    peak_2_peak[[ipeak]] <- use_peaks

    peak_mz <- extract_single(peak_data, use_var = observed_mz, sample_peak = sample_peak, sample = sample)
    mz_matrix[ipeak, names(peak_mz)] <- peak_mz

    peak_height <- extract_single(peak_data, use_var = height,
                                  sample_peak = sample_peak, sample = sample)
    height_matrix[ipeak, names(peak_height)] <- peak_height

    peak_imf <- extract_multiple(peak_data, use_var = imf, sample = sample)
    imf_matrix[ipeak, names(peak_imf)] <- peak_imf
    peak_emf <- extract_multiple(peak_data, use_var = emf, sample = sample)
    emf_matrix[ipeak, names(peak_emf)] <- peak_emf
  }

  n_imf <- purrr::map_int(peak_2_imf, length)
  all_imf <- data.frame(peak = rep(names(peak_2_imf), n_imf), imf = unlist(peak_2_imf),
                        stringsAsFactors = FALSE)
  names(all_imf)[2] <- imf
  all_imf <- dplyr::left_join(all_imf, unique(all_assignments[, c(imf, emf, other_cols)]), by = imf)

  return(list(mz = mz_matrix, height = height_matrix,
              imf = imf_matrix, emf = emf_matrix,
              peak_info = all_imf,
              tic = get_tic(assigned_data)))
}

extract_single <- function(assignment_data, use_var = NULL, sample_peak = "sample_peak", sample = "sample"){

  if (is.null(use_var)) {
    stop("The variable to extract (use_var) must be defined!")
  }

  tmp_data <- unique(assignment_data[, c(use_var, sample_peak, sample)])

  if (nrow(tmp_data) == length(unique(tmp_data[, sample]))) {
    out_data <- tmp_data[, use_var]
    names(out_data) <- tmp_data[, sample]
  } else {
    warning("Data doesn't match samples!")
    out_data <- NA
  }

  out_data

}

extract_multiple <- function(assignment_data, use_var = NULL, sample = "sample"){

  if (is.null(use_var)) {
    stop("The variable to extract (use_var) must be defined!")
  }

  tmp_data <- split(assignment_data[, use_var], assignment_data[, sample])
  tmp_data
}

#' choose a peak
#'
#' When there are multiple peaks in a sample, we need to choose one based on
#' the chosen representative IMF. See **Details** for more information on how
#' that decision is made.
#'
#' @param assignment_data the data.frame of assignment data
#' @param keep_imf which IMF was chosen from the assignments
#' @param imf_mz the assigned IMF M/Z
#' @param imf which variable are the IMFs listed
#' @param sample_peak which variable contains the sample peak information
#' @param sample which variable has the sample information
#'
#' @details Given the `assignment_data` and which IMF to keep, picks which
#' peak should be kept for a given sample. This is made in order of:
#' 1. Does one of the peaks have a **Primary** assignment matching the IMF?
#' 1. Does one of the peaks have a **Secondary** assignment matching the IMF?
#' 1. Which peak has the **smallest** M/Z difference to the chosen IMF?
#'
#' @export
#'
#' @return data.frame of assignment data, with a single peak in each sample.
choose_single_peak <- function(assignment_data, keep_imf, imf_mz = NULL, imf = "IMF",
                               sample_peak = "sample_peak", sample = "sample"){
  split_sample <- split(assignment_data, assignment_data[, sample])

  n_peak <- purrr::map_int(split_sample, function(x){length(unique(x[, sample_peak]))})

  if (all(n_peak == 1)) {
    return(assignment_data)
  }

  split_multi <- split_sample[n_peak > 1]

  decide_peak <- function(sample_data, keep_imf, imf_mz = NULL, imf = "IMF",
                          sample_peak = "sample_peak", sample = "sample"){
    match_imf <- sample_data[, imf] %in% keep_imf

    sample_peak_match_imf <- unique(sample_data[match_imf, sample_peak])

    # assuming there were matches to the chosen IMF
    # return the single matching peak (hopefully)
    if ((length(sample_peak_match_imf) == 1) && (all(!is.na(sample_peak_match_imf)))) {
      return(sample_data[(sample_data[, sample_peak] %in% sample_peak_match_imf), ])


    } else if ((length(sample_peak_match_imf) > 1) && (all(!is.na(sample_peak_match_imf)))) {
      use_loc <- (sample_data[, sample_peak] %in% sample_peak_match_imf) & match_imf
      tmp_data <- sample_data[use_loc, ]
      type_matches <- tmp_data[, "Type"]

      # return the one that matched with a "Primary" assignment
      if ((length(unique(type_matches)) == 2) && (all(!is.na(type_matches)))) {
        use_peak <- tmp_data[type_matches %in% "Primary", imf]
        return(sample_data[(sample_data[, sample_peak] %in% use_peak), ])
      } else if ((length(unique(type_matches)) == 1) && (all(!is.na(type_matches)))) {
        # return the one that had the smallest difference in mass
        mz_diff <- abs(tmp_data$ObservedMZ.Mean - tmp_data$Assigned.M.Z)
        use_peak <- tmp_data[which.min(mz_diff), sample_peak]
        return(sample_data[(sample_data[, sample_peak] %in% use_peak), ])
      }

    }

    # No IMF matches, now what??
    # do it based on closest overall m/z diff to the chosen IMF, which was also
    # passed in, just in case
    mz_diff <- purrr::map_df(split(sample_data, sample_data[, sample_peak]), function(in_peak){
      tmp_frame <- data.frame(mz_diff = abs(in_peak$ObservedMZ.Mean[1] - imf_mz))
      tmp_frame[[2]] <- in_peak[1, sample_peak]
      names(tmp_frame)[2] <- sample_peak
      tmp_frame
    })
    use_peak <- mz_diff[which.min(mz_diff$mz_diff), sample_peak]

    return(sample_data[(sample_data[, sample_peak] %in% use_peak), ])
  }

  multi_data <- purrr::map(split_multi, decide_peak, keep_imf, imf_mz,
                           imf = imf, sample_peak = sample_peak, sample = sample)
  all_data <- purrr::map_df(c(multi_data, split_sample[n_peak == 1]), function(x){x})

  all_data
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
choose_imfs <- function(assignment_data, sample_peak = "sample_peak", imf = "IMF",
sample = "sample", e_value = "IMF.E.Value", min_e_value = 1e-18){

  # easy cases, there is only 1 IMF
  if (nrow(assignment_data) == 1) {
    return(assignment_data[, imf])
  }

  if (length(unique(assignment_data[, imf])) == 1) {
    return(assignment_data[, imf][1])
  }

  # then there are multiple IMFs to choose from
  assignment_single_e_value <- pick_single_evalue(assignment_data,
                                                  sample_peak = sample_peak,
                                                  imf = imf,
                                                  sample = sample,
                                                  e_value = e_value)

  split_assignment_data <- split(assignment_single_e_value, assignment_single_e_value[, imf])
  n_sample <- length(unique(assignment_single_e_value[, sample]))

  imf_stats <- purrr::map_df(split_assignment_data, function(in_imf){
    data.frame(imf = unique(in_imf[, imf]),
               ratio = nrow(in_imf) / n_sample,
               evalue = median_e_values(in_imf[, e_value], min_e_value),
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

  min_evalue <- min(imf_stats$evalue)
  evalue_diffs <- abs(log10(min_evalue) - log10(imf_stats$evalue))
  out_evalue <- !(evalue_diffs > 3)

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

pick_single_evalue <- function(assignment_data, sample_peak = "sample_peak", imf = "IMF",
                               sample = "sample", e_value = "IMF.E.Value"){
  split_imf <- split(assignment_data, paste0(assignment_data[, imf], assignment_data[, sample]))

  assignment_data_single <- purrr::map_df(split_imf, function(x){
    nrow_unique <- nrow(unique(x[, c(imf, e_value)]))
    if (nrow_unique > 1) {
      return(x[which.min(x[, e_value]), ])
    } else {
      return(x[1, ])
    }
  })
}


median_e_values <- function(e_values, min_e_value = 1e-18){
  # handle na
  n_value <- length(e_values)
  if (sum(is.na(e_values)) == n_value) {
    return(1)
  } else {
    e_values <- e_values[!is.na(e_values)]
  }
  e_values[e_values == 0] <- min_e_value
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
create_sudo_peaks <- function(in_assignments, sample_peak = "sample_peak", imf = "IMF"){
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
    tmp_imf <- in_assignments[in_assignments[, sample_peak] %in% tmp_peak, imf]

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

    if (is.null(single_peak$Peak)) {
      peak_id <- ipeak
    } else {
      peak_id <- single_peak$Peak
    }
    if (is.na(single_peak$Sample) & !is.null(single_peak$Scan)) {
      scan <- as.character(single_peak$Scan)
    } else {
      scan <- single_peak$Sample
    }
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
  assign_types <- purrr::map_chr(in_assignment, "Type")

  if (!include_non_primary) {
    in_assignment <- in_assignment[assign_types %in% "Primary"]
  }

  in_assignment <- purrr::map(in_assignment, replace_null)

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
