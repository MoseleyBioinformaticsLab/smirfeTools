#' read SMIRFE assignments
#'
#' Given a SMIRFE JSON output file, read it in, along with other useful
#' information.
#'
#' @param smirfe_assignment the set of assignment results
#' @param .pb a progress bar object
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

#' create pseudo peak
#'
#' Given a master table of assignments, creates a set of "pseudo" peaks
#' that each contain the set of `IMF`s and the `sample_peak` identifiers
#' to pull out of a master table of peak assignments.
#'
#' @param all_assignments the table of all assignments from all samples
#' @param sample_peak which variable holds the sample peak
#' @param imf which variable holds the IMF information
#'
#' @export
#' @return list of pseudo peaks
create_sudo_peaks <- function(all_assignments, sample_peak = "sample_peak", imf = "IMF"){
  all_assignments$grabbed <- FALSE

  peak_index <- 1

  # maximum amount of peaks, likely much less, but this sets an upper bound for
  # us.
  sudo_peaks <- vector("list", length(unique(all_assignments[!all_assignments$grabbed, "sample_peak"])))

  # while there is something grab, go in and get the first peak. We then grab
  # the IMFs for that sample peak, then all the sample peaks for those IMFs,
  # continuing until the set of IMFs and sample peaks no longer change.
  #
  # This should be equivalent to doing set intersections + union across a list
  # of IMFs for sample peaks, but this is much, much faster than checking all
  # of the sample peaks
  while (sum(!all_assignments$grabbed) > 0) {
    tmp_peak <- unique(all_assignments[!all_assignments$grabbed, sample_peak])[1]
    tmp_imf <- all_assignments[all_assignments[, sample_peak] %in% tmp_peak, imf]

    all_peak_imf <- unique(all_assignments[all_assignments[, imf] %in% tmp_imf, sample_peak])
    all_imf_peak <- unique(all_assignments[all_assignments[, sample_peak] %in% all_peak_imf, imf])

    match_peaks <- FALSE

    n_iter <- 1
    while (!match_peaks) {
      all_peak_imf_2 <- unique(all_assignments[all_assignments[, imf] %in% all_imf_peak, sample_peak])
      all_imf_peak_2 <- unique(all_assignments[all_assignments[, sample_peak] %in% all_peak_imf_2, imf])

      match_peaks <- isTRUE(all.equal(all_peak_imf_2, all_peak_imf)) &&
        isTRUE(all.equal(all_imf_peak_2, all_imf_peak))

      all_peak_imf <- all_peak_imf_2
      all_imf_peak <- all_imf_peak_2
      n_iter <- n_iter + 1
    }
    #print(n_iter)
    sudo_peaks[[peak_index]] <- list(imfs = unique(all_imf_peak), peaks = unique(all_peak_imf))

    all_assignments[all_assignments[, imf] %in% all_imf_peak, "grabbed"] <- TRUE
    peak_index <- peak_index + 1

  }
  null_sudos <- purrr::map_lgl(sudo_peaks, is.null)
  sudo_peaks <- sudo_peaks[!null_sudos]
  sudo_peaks
}


#' filter assignments
#'
#' Given an assignment, filter them based on particular criteria, like removing
#' all of the secondary assignments, for example.
#'
#' @param assignment_table the data.frame of assignments
#' @param remove_secondary_assignments should secondary assignments be removed?
#'
#' @return filtered data.frame
#'
#' @export
filter_assignments <- function(assignment_table, remove_secondary_assignments = TRUE){
  if (remove_secondary_assignments) {
    assignment_table <- assignment_table[(assignment_table$Type %in% "Primary"), ]
  }
  assignment_table
}

#' choose a single peak
#'
#' Need to choose a single peak to associate with each IMF. This tries to make
#' a choice either by assignment type, IMF E value, or M/Z difference, in that
#' order.
#'
#' @param peak_frame the data.frame with the IMF information
#'
#' @export
#' @return data.frame with only a single peak
#'
choose_peak <- function(peak_frame){
  primary_frame <- peak_frame[peak_frame$Type %in% "Primary", ] # in case it is due to having secondary assignments
  if (nrow(primary_frame) > 0) {
    uniq_frame <- primary_frame[!duplicated(primary_frame[, c("peak_id", "IMF")]), ]
  } else {
    uniq_frame <- peak_frame[!duplicated(peak_frame[, c("peak_id", "IMF")]), ]
    if (all(is.na(uniq_frame$IMF.E.Value))) {
      uniq_frame$IMF.E.Value <- 1
    }
  }

  if (length(uniq_frame$peak_id) == 1) {
    peak_frame <- peak_frame[peak_frame$peak_id %in% uniq_frame$peak_id, ]
  } else if (sum(uniq_frame$IMF.E.Value == max(uniq_frame$IMF.E.Value)) > 1) {
    mz_diff <- abs(uniq_frame$ObservedMZ.Mean - uniq_frame$Assigned.M.Z)
    peak_frame <- peak_frame[peak_frame$peak_id %in% uniq_frame$peak_id[which.min(mz_diff)], ]
  } else {
    peak_frame <- peak_frame[peak_frame$peak_id %in% uniq_frame$peak_id[which.max(uniq_frame$IMF.E.Value)], ]
  }
  peak_frame
}

count_peaks <- function(x){
  length(unique(x$peak_id))
}

#' single peaks
#'
#' given a data.frame of peak assignments, go through and choose a single peak
#' where needed. Trying to avoid running \code{choose_peak} on every single IMF.
#'
#' @param peak_assignments data.frame of peak assignments
#'
#' @return list of IMFs, each with a single peak
#' @export
one_peak_from_imfs <- function(peak_assignments){
  split_imfs <- split(peak_assignments, peak_assignments$IMF)
  n_peak <- purrr::map_dbl(split_imfs, count_peaks)

  single_peaks <- split_imfs[n_peak == 1]
  other_peaks <- purrr::map(split_imfs[n_peak > 1], choose_peak)

  all_peaks <- c(single_peaks, other_peaks)
  all_peaks
}

peaks_2_imfs <- function(assignment_table){
  split(assignment_table, assignment_table$peak_id)
}

imfs_2_peaks <- function(assignment_table){
  split(assignment_table, assignment_table$IMF)
}

count_peaks <- function(assignment_table){
  length(unique)
}

find_within_sample_peaks <- function(assignment_table){
  split_imfs <- imfs_2_peaks(assignment_table)

  n_peak <- purrr::map_int(split_imfs, function(x){length(unique(x$peak_id))})

  if (sum(n_peak > 1) == 0) {
    out_peaks <- peaks_2_imfs(assignment_table)
    names(out_peaks) <- NULL
    return(out_peaks)
  } else {
    # do something fancy here
    warning("multiple peaks / IMF observed")
  }

}

#' within id
#'
#' for a list of samples, generate a unique group id within each sample
#'
#' @param sample_data list of sample data.frames
#' @param uniq_id which column that defines a unique peak
#' @param merge_id which column to grab for merging on later
#'
#' @export
#' @return list of new data.frames
within_id <- function(sample_data, uniq_id = "peak_id", merge_id = "IMF"){
  lapply(sample_data, function(in_sample){
    grab_id <- unique(in_sample[, uniq_id])
    out_grab <- seq_along(grab_id)
    names(out_grab) <- grab_id
    tmp_data <- in_sample[, c(uniq_id, merge_id)]
    names(tmp_data) <- c("uniq_id", "merge_id")
    tmp_data$group_id <- 0
    tmp_data$sample_loc <- 0
    for (in_id in grab_id){
      all_loc <- tmp_data[, "uniq_id"] %in% in_id
      tmp_data[all_loc, "group_id"] <- out_grab[as.character(in_id)]
      tmp_data[all_loc, "sample_loc"] <- which(all_loc)
    }
    tmp_data
  })
}

#' generate master goups
#'
#' given the groups generated by \code{within_id}, create a master list of groups
#' (peaks) that are present across samples
#'
#' @param within_groups list of data.frames from \code{within_id}
#' @param progress tell the user that data is being processed (default = TRUE)
#'
#' @export
#' @return data.frame of master groups
generate_master <- function(within_groups, .pb = NULL){
  master_groups <- within_groups[[1]]

  not_equal <- TRUE
  round_count <- 1
  init_group <- master_groups

  while (not_equal) {
    if (!is.null(.pb)) {
      tmp_pb <- .pb$clone(deep = TRUE)
    } else {
      tmp_pb <- NULL
    }
    for (i_sample in seq(1, length(within_groups))) {

      knitrProgressBar::update_progress(tmp_pb)


      split_merge <- split(within_groups[[i_sample]][, "merge_id"], within_groups[[i_sample]][, "group_id"])

      for (i_split in seq(1, length(split_merge))) {
        #print(paste("split: ", i_split, sep = ""))
        tmp_split <- split_merge[[i_split]]
        loc_match <- match(tmp_split, master_groups[, "merge_id"])

        na_index <- which(is.na(loc_match))
        na_id <- tmp_split[na_index]
        loc_match <- loc_match[-na_index]

        # Three scenarios:
        # 1: none matching, need to add a new group, incrementing the group id
        # 2: some matching, some not, need to add the ones that were not previously there
        # to the same group
        # 3: some matching to two groups, need to collapse the two in master group into
        # one group

        if ((length(loc_match) == 0) & (length(na_index) > 0)) {

          new_grp_id <- max(master_groups$group_id) + 1
          master_groups <- rbind(master_groups,
                                 data.frame(uniq_id = NA, merge_id = na_id,
                                            group_id = new_grp_id, sample_loc = NA,
                                            stringsAsFactors = FALSE))

        } else if ((length(loc_match) > 0) & (length(na_index) > 0)) {
          master_grp_id <- master_groups[loc_match, "group_id"]
          uniq_grp_id <- unique(master_grp_id)

          if (length(uniq_grp_id) > 1) {
            #stop("merging groups")
            master_loc <- master_groups$group_id %in% uniq_grp_id
            new_grp_id <- max(master_groups$group_id) + 1
            master_groups$group_id[master_loc] <- new_grp_id
            #message(cat("Merging groups", uniq_grp_id, sep =" "))

          } else {
            new_grp_id <- uniq_grp_id
          }

          master_groups <- rbind(master_groups,
                                 data.frame(uniq_id = NA, merge_id = na_id,
                                            group_id = new_grp_id, sample_loc = NA,
                                            stringsAsFactors = FALSE))
        }
      }
    }
    # check if what we have now is the same as before
    # if there are no changes, then we can stop
    # otherwise we loop through the samples again, and see
    # if anything changes
    if (identical(master_groups, init_group)) {
      not_equal <- FALSE
    } else {
      init_group <- master_groups
      round_count <- round_count + 1
    }
  }

  # if (progress) {
  #   message("Done!")
  # }

  # finally, split on the group id, and then renumber them and do a proper make.names
  # turns out we don't use the uniq_id or sample_loc anywhere later, so we can
  # discard them. All matching in samples is done using the merge_id
  split_groups <- split(master_groups$merge_id, master_groups$group_id)
  names(split_groups) <- make.names(seq_len(length(split_groups)))

  split_groups
}


#' extract peaks
#'
#' To be computably useful, we need to associate each  the peak information for each IMF needs to be in a
#' matrix type format. This takes a list of peak-IMF data.frames, and generates
#' a list of matrices to hold the data, where the matrix has IMF rows and sample
#' columns.
#'
#' @param assignment_list list of assignments, one for each sample, generated by \code{read_smirfe_assignment}
#' @param .pb a progress bar object from either dplyr or knitrProgressBar
#'
#' @return list with height, emf, peak, and type, and emf_imf list
#' @export
#'
extract_peaks <- function(assignment_list, .pb = NULL){
  get_imfs <- function(x){x$assignments$IMF}
  get_tic <- function(x){x$tic}
  get_sample <- function(x){x$sample}
  names(assignment_list) <- purrr::map_chr(assignment_list, get_sample)

  filtered_assignments <- purrr::map(assignment_list, function(x){
    filter_assignments(x$assignments)
  })

  n_peaks <- purrr::map_int(filtered_assignments, function(x){
    length(unique(x$peak_id))
  })
  total_peaks <- round(sum(n_peaks) * .9)

  peak_list <- vector("list", total_peaks)

  possible_peaks <- find_within_sample_peaks(filtered_assignments[[1]])

  peak_imf_map <- purrr::map(possible_peaks, function(x){
    x$IMF
  })

  for (isample in filtered_assignments) {

  }

  all_imfs <- purrr::map(assignment_list, get_imfs)
  all_imfs <- unique(unlist(all_imfs))
  imf_index <- data.frame(IMF = all_imfs, stringsAsFactors = FALSE)

  height_matrix <- matrix(NA, nrow = nrow(imf_index), ncol = length(assignment_list))
  emf_matrix <- matrix(list(), nrow = nrow(imf_index), ncol = length(assignment_list))
  peak_matrix <- matrix(NA, nrow = nrow(imf_index), ncol = length(assignment_list))
  type_matrix <- matrix("NA", nrow = nrow(imf_index), ncol = length(assignment_list))
  emf_imf_mappings <- vector(mode = "list", length = length(assignment_list))
  names(emf_imf_mappings) <- names(assignment_list)

  rownames(height_matrix) <- rownames(emf_matrix) <- rownames(peak_matrix) <- rownames(type_matrix) <- all_imfs
  colnames(height_matrix) <- colnames(emf_matrix) <- colnames(peak_matrix) <- colnames(type_matrix) <- names(assignment_list)

  for (isample in names(assignment_list)) {
    if (!is.null(.pb)) {
      knitrProgressBar::update_progress(.pb)
    }

    imfs <- one_peak_from_imfs(assignment_list[[isample]]$assignments)

    for (iimf in seq_along(imfs)) {
      #print(iimf)
      use_imf <- imfs[[iimf]]
      height_matrix[use_imf$IMF[1], isample] <- use_imf$Height.Mean[1]
      emf_matrix[use_imf$IMF[1], isample][[1]] <- use_imf$EMF
      peak_matrix[use_imf$IMF[1], isample] <- use_imf$peak_id[1]
      type_matrix[use_imf$IMF[1], isample] <- use_imf$Type[1]
    }

    joined_imfs <- purrr::map_df(imfs, function(x){x})
    emf_imf_mappings[[isample]] <- split(joined_imfs$IMF, joined_imfs$EMF)
  }
  list(height = height_matrix,
       emf = emf_matrix,
       peaks = peak_matrix,
       type = type_matrix,
       tic = purrr::map_dbl(assignment_list, get_tic),
       samples = names(assignment_list),
       emf_imf = emf_imf_mappings)
}


unique_nona <- function(value){unique(value[!is.na(value)])}

test_na <- function(in_vector, n = 1){
  length(unique_nona(in_vector)) <= n
}

test_na_matrix <- function(in_matrix){
  unique_na <- purrr::map_lgl(seq(1, ncol(in_matrix)), function(icol){
    test_na(in_matrix[, icol])
  })
  all(unique_na)
}

remove_na_cols <- function(in_matrix){
  not_na <- purrr::map_lgl(seq(1, ncol(in_matrix)), function(icol){
    !all(is.na(in_matrix[, icol]))
  })
  in_matrix[, not_na, drop = FALSE]
}

clean_imfpeak_matrix <- function(imf_2_peaks, n = 1){
  # find any sample that actually has more than one peak in it
  unique_na <- purrr::map_lgl(seq(1, ncol(imf_2_peaks)), function(icol){
    test_na(imf_2_peaks[, icol], n = n)
  })
  # remove it
  imf_2_peaks <- imf_2_peaks[, unique_na, drop = FALSE]

  # after removing a sample with multiple peaks to our IMFs, then make sure we
  # haven't made any of our rows all NA, not likely, but could happen
  na_rows <- purrr::map_lgl(seq(1, nrow(imf_2_peaks)), function(irow){
    all(is.na(imf_2_peaks[irow, ]))
  })

  imf_2_peaks <- imf_2_peaks[!na_rows, , drop = FALSE]
  imf_2_peaks
}

filter_extracted <- function(keep_imfs, extracted_data, imf_2_imf){
  # these ones are simple to extract
  filter_vars <- c("height", "peaks", "type")
  extracted_data[filter_vars] <- purrr::map(extracted_data[filter_vars], function(in_var){
    in_var[keep_imfs, ]
  })

  # but emf needs some more work
  names(imf_2_imf) <- keep_imfs

  emfs <- extracted_data$emf
  emfs2 <- emfs[keep_imfs, ]

  for (isample in colnames(emfs)) {
    for (imf in keep_imfs) {
      out_emfs <- unique(unlist(emfs[imf_2_imf[[imf]], isample]))
      if (!is.null(out_emfs)) {
        emfs2[imf, isample][[1]] <- out_emfs
      }

    }
  }

  extracted_data$emf <- emfs2
  extracted_data
}

#' find duplicate peaks
#'
#' @param extracted_data generated by \code{extract_peaks}
#' @param .pb a progress object, should be equal to the number of samples used.
#'
#' @export
#' @return deduplicated list
find_duplicate_peaks <- function(extracted_data, .pb = NULL){
  split_imfs <- purrr::map(seq(1, ncol(extracted_data$peaks)), function(in_sample){
    split(rownames(extracted_data$peaks), extracted_data$peaks[, in_sample])
  })
  names(split_imfs) <- colnames(extracted_data$peaks)

  all_imfs <- data.frame(imf = rownames(extracted_data$height), grabbed = FALSE,
                         stringsAsFactors = FALSE)
  all_imfs$others <- vector("list", length = nrow(all_imfs))

  for (isplit in split_imfs) {
    if (!is.null(.pb)) {
      knitrProgressBar::update_progress(.pb)
    }
    isplit_diffs <- purrr::map(isplit, function(x){
      base::setdiff(x, all_imfs$imf[all_imfs$grabbed])
    })
    isplit_null <- purrr::map_lgl(isplit_diffs, function(x){length(x) == 0})
    isplit_diffs <- isplit_diffs[!isplit_null]

    for (ipeak in names(isplit_diffs)) {
      all_peaks <- extracted_data$peaks[isplit_diffs[[ipeak]], , drop = FALSE]
      all_peaks <- remove_na_cols(all_peaks)

      curr_imfs <- purrr::map(colnames(all_peaks), function(isample){
        use_peak <- as.character(unique_nona(all_peaks[, isample]))
        split_imfs[[isample]][use_peak]
      })
      curr_imfs <- unique(unlist(curr_imfs))

      imf_2_peaks <- extracted_data$peaks[curr_imfs, colnames(all_peaks), drop = FALSE]
      imf_2_peaks <- clean_imfpeak_matrix(imf_2_peaks, n = 2)

      together_imfs <- rownames(imf_2_peaks)

      all_imfs[all_imfs$imf %in% together_imfs, "grabbed"] <- TRUE
      all_imfs[all_imfs$imf %in% together_imfs[1], "others"][[1]] <- list(together_imfs)
    }
  }
  all_imfs <- all_imfs[all_imfs$grabbed, ]
  null_imfs <- purrr::map_lgl(all_imfs$others, is.null)
  all_imfs <- all_imfs[!null_imfs, ]

  extracted_data2 <- filter_extracted(all_imfs$imf, extracted_data, all_imfs$others)
  extracted_data2$imf_2_imf <- all_imfs$others
  extracted_data2
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
