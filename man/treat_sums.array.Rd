% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recruit.R
\name{treat_sums.array}
\alias{treat_sums.array}
\title{Sum accrual array by experimental arm (including control).}
\usage{
\method{treat_sums}{array}(
  x,
  control_total = FALSE,
  no_treat_arms,
  shared_control = TRUE,
  na.rm = TRUE
)
}
\arguments{
\item{x}{Accrual array with dimensions Weeks, Arms and Centres.}

\item{control_total}{Logical; if TRUE return single total for all
control arms (not used if \code{shared_control} is TRUE); defaults to FALSE.}

\item{no_treat_arms}{Number of treatment arms (as opposed to control
arms).}

\item{shared_control}{TRUE if all treatment arms share the
same control arm; FALSE if each treatment arm has its own
control. Defaults to TRUE.}
}
\value{
vector of total accrual by experimental arm.
}
\description{
Sum accrual array by experimental arm (including control).
}
