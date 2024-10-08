% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-graphics.R
\name{accrual_plot_from_file}
\alias{accrual_plot_from_file}
\title{Plot predicted recruitment from file containing a CSV from
a single run}
\usage{
accrual_plot_from_file(
  file_prefix = "accrual",
  plot_prefix = "accrual-from-file",
  run_time = "2024-08-07-18-35-09",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/")
)
}
\arguments{
\item{file_prefix}{Consistent beginning of filenames holding
arm closure data. Defaults to \code{accrual}.}

\item{plot_prefix}{Prefix for file name to identify plot type.
Defaults to \code{accrual_plot}.}

\item{run_time}{Specify a particular instance of \code{biomkrAccrual()}
execution using a date-time format \code{yyyy-mm-dd-hh-mm-ss}.
Used to select which files will be summarised.}

\item{output_path}{Directory where the output files from the
\code{biomkrAccrual()} instance are located.}

\item{figs_path}{Folder where figures generated during execution
will be stored; defaults to the \code{figures} subdirectory in
\code{output_path}.}
}
\description{
Plot predicted recruitment from file containing a CSV from
a single run
}
