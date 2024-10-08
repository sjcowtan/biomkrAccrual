% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recruit.R
\name{apply_site_cap}
\alias{apply_site_cap}
\title{Implement site cap on a week's accrual.}
\usage{
apply_site_cap(obj, ...)
}
\arguments{
\item{obj}{Accrual object}

\item{site_cap}{Maximum number of patients per site}
}
\value{
Modified accrual object with capped week's accrual and
with any capped sites removed from active_sites
}
\description{
Implement site cap on a week's accrual.
}
