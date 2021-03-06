% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reading_assignments.R
\name{remove_only_labeled_emfs}
\alias{remove_only_labeled_emfs}
\title{Remove Only Labeled EMFs}
\usage{
remove_only_labeled_emfs(assignment_data, remove_s = TRUE)
}
\arguments{
\item{assignment_data}{a data.frame of assignments for peaks}

\item{remove_s}{logical, should EMFs with Sulfur also be removed?}
}
\value{
data.frame
}
\description{
Removes EMFs that \emph{only} contain labeled IMFs.
}
\details{
When working with metabolomics data where labeled precursors are used
and there is a reasonable expectation that label is incorporated,
it is highly unlikely that there are cases where the un-labeled
IMF among a set of EMFs will not be observed.
However, this is not baked into the
SMIRFE assignment algorithm, so we provide this functionality here.

Note that this should be run \strong{after} \code{read_smirfe_assignment} and
\strong{before} \code{extract_assigned_data}.
}
