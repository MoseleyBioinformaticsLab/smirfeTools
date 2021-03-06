% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reading_assignments.R
\name{create_sudo_peaks}
\alias{create_sudo_peaks}
\title{create pseudo peak}
\usage{
create_sudo_peaks(
  in_assignments,
  sample_peak = "Sample_Peak",
  imf = "complete_IMF"
)
}
\arguments{
\item{in_assignments}{the table of all assignments from all samples}

\item{sample_peak}{which variable holds the sample peak}

\item{imf}{which variable holds the IMF information}
}
\value{
list of pseudo peaks
}
\description{
Given a master table of assignments, creates a set of "pseudo" peaks
that each contain the set of \code{IMF}s and the \code{sample_peak} identifiers
to pull out of a master table of peak assignments.
}
