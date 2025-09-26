#' Prints synthesised recruitment as a dataframe with
#' one column for each treatment and control arm, and
#' one row for each week.
#' 
#' @name print
#' @aliases print.accrual
#' 
#' @param x An object of class `accrual`.
#' @param ... Additional arguments passed to print().
#' 
#' @importFrom S7 new_generic method
#' 
#' @export
#' 
S7::new_generic("print", "accrual")
S7::method(print, accrual) <- function(x, ...) {
  print(data.frame(rowSums(x@accrual, dims = 2)), ...)
}


#' Summary of predicted accrual.
#' 
#' @name summary
#' @aliases summary.accrual
#' 
#' @param x An object of class `accrual`.
#' @param ... Additional arguments passed to summary().
#' 
#' @importFrom S7 new_generic method
#' 
#' @export
#' 
S7::new_generic("summary", "accrual")
S7::method(summary, accrual) <- function(x, ...) {
  
  # Summary of accrual by arm
  cat("Recruitment by experimental arm\n")
  print(summary(data.frame(rowSums(x@accrual, dims = 2)), ...))

  # Summary of accrual by site
  cat("\nRecruitment by site\n")
  print(summary(data.frame(rowSums(
    aperm(x@accrual, c(1, 3, 2)),
    dims = 2
  )), ...))

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
#' @param run_time Specify a particular instance of `biomkrAccrual()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`.
#' @param output_path = Directory where the output files from the 
#' `biomkrAccrual()` instance are located.
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
    figs_path = figs_path,
    target_arm_size = accrual_obj@target_arm_size,
    target_control = accrual_obj@target_control,
    target_interim = accrual_obj@target_interim,
    accrual_period = accrual_obj@accrual_period,
    interim_period = accrual_obj@interim_period
  )
}


#' Prints initial trial structure as a matrix of prevalences by recruitment
#' and experimental arms.
#' 
#' @name print
#' @aliases print.trial_strucutre
#' 
#' @param x An object of class `trial_structure`.
#' @param ... Additional arguments passed to print().
#' 
#' @importFrom S7 new_generic method
#' @importFrom withr with_options
#' 
#' @export
#' 
S7::new_generic("print", "trial_structure")
S7::method(print, trial_structure) <- function(x, ...) {
  
  orig_struct_df <- data.frame(
    x@treatment_arm_struct_start
  )

  colnames(orig_struct_df) <- names(x@treatment_arm_ids_start)
  rownames(orig_struct_df) <- x@recruit_arm_names

  print(orig_struct_df, ...)

}


#' Prints initial trial structure as a matrix of prevalences by recruitment
#' and experimental arms.
#' 
#' @name plot
#' @aliases plot.trial_structure
#' 
#' @param x An object of class `trial_structure`.
#' 
#' @importFrom S7 new_generic method
#' @importFrom stats reshape
#' @importFrom grDevices palette.colors
#' 
#' @export
#' 
S7::new_generic("plot", "trial_structure")
S7::method(plot, trial_structure) <- function(x) {
  
  orig_struct_df <- data.frame(
    x@treatment_arm_struct_start
  )

  colnames(orig_struct_df) <- names(x@treatment_arm_ids_start)
  orig_struct_df$Recruitment <- x@recruit_arm_names

  orig_struct_df <- stats::reshape(
    orig_struct_df,
    direction = "long",
    v.names = "Recruits",
    varying = list(names(x@treatment_arm_ids_start)),
    idvar = "Recruitment",
    timevar = "Treatment",
    times = names(x@treatment_arm_ids_start)
  )

  orig_struct_df$Recruits <- as.factor(as.integer(orig_struct_df$Recruits))

  p <- ggplot2::ggplot(
    data = orig_struct_df,
    ggplot2::aes(x = Treatment, y = Recruitment, fill = Recruits)
  ) +
    ggplot2::geom_tile(
      color = "white",
      lwd = 1.5,
      linetype = 1
    ) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(
      values = c("grey85", grDevices::palette.colors(4)[4]),
      labels = c("No", "Yes")
    ) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(
      x = "Treatment arm",
      y = "Biomarker",
      title = "Trial structure",
      subtitle = paste0(
        ifelse(x@shared_control, "Shared", "Individual"), 
        " control arm",
        ifelse(x@shared_control, "", "s")
      )
    ) +
    theme_bma(base_size = 16)

  return(p)
}


#' Summary of trial structure.
#' 
#' @name summary
#' @aliases summary.trial_structure
#' 
#' @param x An object of class `trial_structure`.
#' @param ..., Additional arguments passed to print().
#' @param digits A non-null value for digits specifies 
#' the minimum number of significant digits to be 
#' printed in values. The default, NULL, uses 
#' getOption("digits"). (For the interpretation for 
#' complex numbers see signif.) Non-integer values 
#' will be rounded down, and only values greater than 
#' or equal to 1 and no greater than 22 are accepted.
#' 
#' @importFrom S7 new_generic method
#' @importFrom withr with_options
#' 
#' @export
#' 
S7::new_generic("summary", "trial_structure")
S7::method(summary, trial_structure) <- 
  function(
    x,
    ...,
    digits = max(3L, getOption("digits") - 3L)
  ) {

    summary_ls <- vector(mode = "list", length = 1)

    # Site prevalences by recruitment arm

    orig_prev_df <- data.frame(
      x@recruit_arm_prevalence_start,
      row.names = x@recruit_arm_names
    )
    colnames(orig_prev_df) <- paste("Site", seq_len(ncol(orig_prev_df)))

    summary_ls$site_prev <- withr::with_options(
      list(scipen = 10),
      print(round(orig_prev_df, digits = digits), ...)
    )

    summary_ls
  }