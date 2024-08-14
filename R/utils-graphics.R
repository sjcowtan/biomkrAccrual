 
#' Theme for ggplot2
#' @param base_size Legend title size, all other sizes scaled appropriately 
#' to this
#' @param base_family Font family
#' 
#' @import ggplot2
#' 
theme_bma <- function(base_size = 10, base_family = "") {

  `%+replace%` <- ggplot2::`%+replace%`

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(size = base_size + 4),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = base_size - 1)
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(l = 0, r = base_size + 1), 
        angle = 90
      ),
      legend.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size + 2),
      strip.background = ggplot2::element_rect(fill = "grey90")
    )
}


#' Generate arm closure summaries
#' 
#' @param file_prefix Consistent beginning of filenames holding 
#' arm closure data. Defaults to `closures`.
#' @param run_time Specify a particular instance of `spine()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`. 
#' Used to select which files will be summarised.
#' @param output_path Directory where the input files are located
#' and the output files will be written.
#' 
#' @export
#' 
#' @importFrom stats sd
#' @importFrom utils read.csv write.csv
#' 
get_arm_closures <- function(
  file_prefix = "closures",
  run_time = "2024-08-07-18-35-09",
  output_path = "../biomkrAccrual_output_data/"
) {
  # What output files do we have?
  filenames <- list.files(
    output_path, 
    pattern = paste0("^", file_prefix, "-", run_time, ".*.csv"), 
    full.names = TRUE
  )

  # Read files
  closures_ls <- lapply(filenames, read.csv)

  print(filenames)

  # Summarise files
  summ <- lapply(closures_ls, function(d) summary(d[, -1]))
  name_types <- c("1", "2", "mixed", "unbalanced", "multirate", "multimix")
  names(summ) <- c(
    paste0("gamma_rate_closures_", name_types),
    paste0("fixed_rate_closures_", name_types)
  )

  print(summ)

  write.csv(
    summ, 
    paste0(output_data, "arm_closure_summary-", run_time, ".csv")
  )

  # Standard Deviations
  sd_mx <- t(sapply(
    closures_ls, 
    function(a_df) {
      sapply(
        seq_len(ncol(a_df)),
        function(i) stats::sd(a_df[i, ], na.rm = TRUE)
      )
    }
  ))

  sd_mx <- sd_mx[, -1]

  colnames(sd_mx) <- paste0("T", 1:3, "_sd")
  rownames(sd_mx) <- c(
    paste0("gamma_rate_closures_", name_types),
    paste0("fixed_rate_closures_", name_types)
  )

  write.csv(
    as.data.frame(sd_mx), 
    paste0(output_data, "arm_closures_sd", run_time, ".csv")
  )

  print(sd_mx)

  return(list(summ, sd_mx))
}


#' Plot a scatter plot with error bars for time to recruit against prevalence
#' 
#' @param prevalences Vector of prevalences
#' @param e_time Vector of expected times
#' @param v_time Vector of variances
#' 
ggscatterError <- function(prevalences, e_time, v_time) {
  plot_df <- data.frame(
    Prevalences = prevalences,
    Months = e_time,
    ymin = e_time - (1.96 * sqrt(v_time)),
    ymax = e_time + (1.96 * sqrt(v_time))
  )

  plot_col <- grDevices::palette.colors(3)[3]

  ggplot2::ggplot(
    plot_df, 
    ggplot2::aes(x = factor(get("Prevalences")), y = get("Months"))
  ) +
    ggplot2::geom_point(col = plot_col) +
    ggplot2::geom_pointrange(
      col = plot_col,
      ggplot2::aes(ymin = get("ymin"), ymax = get("ymax"))
    ) +
    ggplot2::labs(
      x = "Prevalences"
    ) +
    theme_bma(14)
}


#' Creates sensitivity plot for Poisson distribution with 
#' specified site rates for simultaneous start
#' 
#' @param plot_name Name for sensitivity plot file (.png will be appended).
#' @param target_arm_size Number of patients to be recruited.
#' @param site_rates Total site rate.
#' 
do_sensitivity_plot_simultaneous <- function( 
  target_arm_size, 
  site_rates,
  plot_name,
  figs_path
) {
  # Fixed, regularly spaced prevalences
  prevalences <- seq(0.1, 0.9, by = 0.1)

  # Expectation
  e_time <- target_arm_size / (prevalences * site_rates)

  # Variance
  v_time <- target_arm_size / (prevalences * site_rates)^2

  p <- ggscatterError(prevalences, e_time, v_time)

  print(p)

  ggplot2::ggsave(paste0(figs_path, plot_name, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
  )
}


#' Plot predicted recruitment from file containing a CSV from
#' a single run
#' 
#' @param file_prefix Consistent beginning of filenames holding 
#' arm closure data. Defaults to `accrual`.
#' @param plot_prefix Prefix for file name to identify plot type. 
#' Defaults to `accrual_plot`.
#' @param run_time Specify a particular instance of `spine()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`. 
#' Used to select which files will be summarised.
#' @param output_path Directory where the output files from the 
#' `spine()` instance are located.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' 
#' @export
#' 
accrual_plot_from_file <- function(
  file_prefix = "accrual",
  plot_prefix = "accrual_plot",
  run_time = "2024-08-07-18-35-09",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/")
) {
  # Validate input

  checkmate::assert_directory_exists(
    file.path(output_path), 
    access = "rx"
  )

  input_file <- paste0(
    output_path, file_prefix, "-", run_time, ".csv"
  )

  checkmate::assert_file_exists(file.path(input_file))

  makeifnot_dir(figs_path, min_access = "rwx")

  accrual_raw_df <- utils::read.csv(file.path(input_file))

  # Get unique arm identifiers
  arm_names <- unique(sapply(
    strsplit(names(accrual_raw_df), "\\."), 
    getElement,
    1
  ))

  # Make dataframe of arm recruitment per week by summing site-arm columns
  accrual_df <- as.data.frame(
    lapply(
      arm_names, 
      function(n) rowSums(accrual_raw_df[startsWith(names(accrual_raw_df), n)])
    ),
    col.names = arm_names
  )

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


#' Convert accrual data to long format dataframe of class
#' `accrualplotformat`.
#' 
#' @param accrual_df Wide format accrual data.
#' 
#' @importFrom stats reshape
#' 
accrual_to_long <- function(accrual_df) {
  # Convert to cumulative sums
  accrual_df <- cumsum(accrual_df)

  arm_names <- names(accrual_df)

  # Add week information
  accrual_df$Week <- seq_len(nrow(accrual_df))

  # Convert to long format
  accrual_df <- stats::reshape(
    accrual_df,
    direction = "long",
    varying = list(arm_names),
    timevar = "Arm",
    times = arm_names,
    v.names = "Recruitment",
    idvar = "Week"
  )

  accrual_df$Arm <- factor(accrual_df$Arm)

  # Define an S3 class so we can have a custom plot command
  class(accrual_df) <- c("accrualplotdata", class(accrual_df))

  return(accrual_df)
}


#' S3 method to plot predicted recruitment from a long format 
#' dataframe of class "accrualplotdata".
#' 
#' @param data long format dataframe with columns "Week", 
#' "Arm" and "Recruitment".
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
#' @import ggplot2
#' @importFrom grDevices palette.colors
#' @importFrom rlang .data
#' 
plot.accrualplotdata <- function(
  data,
  plot_prefix = "accrual_plot",
  run_time = NULL,
  figs_path = paste0(output_path, "figures/")
) {
  
  accrual_df <- data
  arm_names <- levels(accrual_df$Arm)

  p <- ggplot2::ggplot(
    accrual_df, 
    ggplot2::aes(
      x = .data$Week, 
      y = .data$Recruitment, 
      group = .data$Arm, 
      color = .data$Arm
    )
  ) +
    ggplot2::geom_line() +
    # Use colourblind friendly Okabe-Ito palette
    ggplot2::scale_colour_manual(
      values = grDevices::palette.colors(length(arm_names))
    ) +
    theme_bma(base_size = 12)

  ggplot2::ggsave(
    paste0(figs_path, plot_prefix, "-", run_time, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
  )

  print(p)
}