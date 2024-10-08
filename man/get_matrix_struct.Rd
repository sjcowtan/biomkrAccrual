% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prevalence.R
\name{get_matrix_struct}
\alias{get_matrix_struct}
\title{Converts trial structure and prevalence information into matrix form}
\usage{
get_matrix_struct(arms_ls, recruit_arm_prevalence)
}
\arguments{
\item{arms_ls}{List of lists of recruitment arms which recruit to
each treatment arm.}

\item{recruit_arm_prevalence}{Matrix of prevalences with one row per
biomarker and one column per site; each column sums to 1.}
}
\value{
Logical matrix with one row per biomarker and one column per
treatment arm; TRUE where the treatment recruits from that biomarker.
}
\description{
Converts trial structure and prevalence information into matrix form
}
