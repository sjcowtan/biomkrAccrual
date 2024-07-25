 
#' Custom colours (colourblind friendly) -
#' omits colours nearest black and white
#' @param palette Defaults to "viridis"
#' @param no_colours Defaults to 4
#' 
#' @importFrom viridisLite viridis
#' 
bma_colours <- function(palette = "viridis", no_colours = 4) {
  viridisLite::viridis(no_colours * 2, option = palette)[
    seq(2, no_colours * 2, by = 2)
  ]
}


#' Theme for ggplot2
#' @param base_size Legend title size, all other sizes scaled appropriately 
#' to this
#' @param base_family Font family
#' 
#' @import ggplot2
#' 
theme_bma <- function(base_size = 10, base_family = "") {
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
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size),
      strip.background = ggplot2::element_rect(fill = "grey90")
    )
}


#' Display arm closure summaries
#' @param file_prefix Consistent beginning of filename
#' @param output_path Directory where output .csvs are written
#' @export
#' @importFrom stats sd
#' @importFrom utils read.csv write.csv
ggclosures <- function(
  file_prefix = "arm_closures_24-04-29_",
  output_path = "output_data/"
) {
  # What output files do we have?
  filenames <- list.files(
    output_path, 
    pattern = paste0("^", file_prefix, ".*.csv"), 
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

  write.csv(summ, "output_data/arm_closure_summary.csv")

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

  write.csv(as.data.frame(sd_mx), "output_data/arm_closures_sd.csv")

  print(sd_mx)
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

  plot_col <- bma_colours()[1]

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