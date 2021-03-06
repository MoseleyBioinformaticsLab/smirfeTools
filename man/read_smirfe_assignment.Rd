% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reading_assignments.R
\name{read_smirfe_assignment}
\alias{read_smirfe_assignment}
\title{read SMIRFE assignments}
\usage{
read_smirfe_assignment(
  smirfe_assignment,
  assigned_only = TRUE,
  sample_id = NULL,
  .pb = NULL
)
}
\arguments{
\item{smirfe_assignment}{the set of assignment results}

\item{assigned_only}{whether to return peaks with assignment only}

\item{sample_id}{a user provided sample id to override what is in the metadata}

\item{.pb}{a progress bar object}
}
\value{
list of tic, assignments, sample
}
\description{
Given a SMIRFE JSON output file, read it in, along with other useful
information.
}
