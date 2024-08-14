#' Prints synthesised recruitment as a dataframe with
#' one column for each treatment and control arm, and
#' one row for each week.
#' 
#' @name print
#' @aliases print.accrual
#' 
#' @param x An object of class `accrual`.
#' 
#' @importFrom S7 new_generic method
#' 
#' @export
#' 
S7::new_generic("print", "accrual")
S7::method(print, accrual) <- function(x) {
  print(data.frame(rowSums(x@accrual, dims = 2)))
}


#' Summary of predicted accrual.
#' 
#' @name summary
#' @aliases summary.accrual
#' 
#' @param x An object of class `accrual`.
#' 
#' @importFrom S7 new_generic method
#' 
#' @export
#' 
S7::new_generic("summary", "accrual")
S7::method(summary, accrual) <- function(x) {
  
  # Summary of accrual by arm
  cat("Recruitment by experimental arm\n")
  print(summary(data.frame(rowSums(x@accrual, dims = 2))))

  # Summary of accrual by site
  cat("\nRecruitment by site\n")
  print(summary(data.frame(rowSums(
    aperm(x@accrual, c(1, 3, 2)),
    dims = 2
  ))))

  # Summary of phase change weeks
  cat("\nExperimental arm closure weeks\n")
  acw <- as.vector(x@phase_changes)
  names(acw) <- dimnames(x@accrual)$Arms[seq_len(length(acw))]
  print(acw)

  # Summary of accrual totals by arm
  cat("\nAccrual totals by experimental arm\n")
  print(treat_sums(x))
}


#' Plot method for an object of class `accrual`.  Creates
#' a line plot of cumulative recruitment, grouped by trial arm,
#' using ggplot2.
#' 
#' @name plot
#' @aliases plot.accrual
#' 
#' @param accrual_obj Object of class `accrual`.
#' @param plot_prefix Prefix for file name to identify plot type.
#' Defaults to `accrual_plot`.
#' @param run_time Specify a particular instance of `spine()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`.
#' @param output_path = Directory where the output files from the 
#' `spine()` instance are located.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' 
#' @export
#' 
S7::new_generic("plot", "accrual")
S7::method(plot, accrual) <- function(
  accrual_obj,
  plot_prefix = "accrual_plot",
  run_time = "2024-08-07-18-35-09",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/")
) {
  accrual_ar <- accrual_obj@accrual

  # Sum across sites
  accrual_df <- data.frame(rowSums(accrual_ar, dims = 2))

  # Convert to long format of class "accrualplotdata"
  accrual_df <- accrual_to_long(accrual_df)

  # Plot and save plot in figs_path
  plot(
    accrual_df, 
    plot_prefix = plot_prefix,
    run_time = run_time,
    figs_path = figs_path
  )
}

