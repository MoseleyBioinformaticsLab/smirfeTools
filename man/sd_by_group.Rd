% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequency_to_mz_cutoff.R
\name{sd_by_group}
\alias{sd_by_group}
\title{sd by group}
\usage{
sd_by_group(grouped_samples, ...)
}
\arguments{
\item{grouped_samples}{list of \code{assigned_data}}

\item{...}{parameters for \code{find_confident_frequency_sd}}
}
\value{
data.frame
}
\description{
Assuming you have some criteria by which you have split samples into a list,
and now you want to calculate the IMF frequency SD for each one independently.
}
