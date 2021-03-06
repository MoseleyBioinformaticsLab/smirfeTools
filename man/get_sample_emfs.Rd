% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emf_voting.R
\name{get_sample_emfs}
\alias{get_sample_emfs}
\title{extract sample level elemental molecular formulas}
\usage{
get_sample_emfs(
  sample_assignments,
  sample_id,
  evalue_cutoff = 0.98,
  use_corroborating = TRUE,
  emf_classifications = NULL
)
}
\arguments{
\item{sample_assignments}{a sample level assignment data.frame from \code{read_smirfe_assignment}}

\item{sample_id}{which sample is it}

\item{evalue_cutoff}{what value should be used to exclude EMFs? (default is 0.98)}

\item{use_corroborating}{should other adducts be added for corroborating evidence (default is \code{TRUE})}

\item{emf_classifications}{data.frame of classifications for isotopologue_EMFs}
}
\value{
list
}
\description{
given the assignment data.frame, group and extract EMF level information
}
