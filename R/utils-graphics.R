 
#' Theme for ggplot2
#' @param base_size Legend title size, all other sizes scaled appropriately 
#' to this
#' 
#' @import ggplot2
#' 
theme_bma <- function(
  base_size = 10
) {

  `%+replace%` <- ggplot2::`%+replace%`

  ggplot2::theme_bw(base_size = base_size) %+replace%
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(
        size = base_size + 6,
        margin = margin(0, 0, 13, 0),
        hjust = 0,
        face = "bold"
      ),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(
        size = base_size + 4,
        margin = margin(0, 0, 13, 0),
        hjust = 0
      ),
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


#' Generate arm closure summaries for batches
#' 
#' @param file_prefix Consistent beginning of filenames holding 
#' arm closure data. Defaults to `closures`.
#' @param run_time Specify a particular instance of `biomkrAccrual()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`. 
#' Used to select which files will be summarised.
#' @param output_path Directory where the input files are located
#' and the output files will be written.
#' @param keep_files Save data files and plots generated during the run. 
#' Defaults to TRUE.
#' 
#' @export
#' 
#' @importFrom stats sd
#' @importFrom utils read.csv write.csv
#' 
##### Nothing is calling this?
get_arm_closures2 <- function(
  file_prefix = "closures",
  run_time = "2024-08-07-18-35-09",
  output_path = "../biomkrAccrual_output_data/",
  keep_files = TRUE
) {
  # What output files do we have?
  filenames <- list.files(
    output_path, 
    pattern = paste0("^", file_prefix, "-", run_time, ".*.csv"), 
    full.names = TRUE
  )

  # Read files
  closures_ls <- lapply(filenames, read.csv)

  # Summarise files
  summ <- lapply(closures_ls, function(d) summary(d[, -1]))
  name_types <- c("1", "2", "mixed", "unbalanced", "multirate", "multimix")
  names(summ) <- c(
    paste0("gamma_rate_closures_", name_types),
    paste0("fixed_rate_closures_", name_types)
  )

  write.csv(
    summ, 
    paste0(output_path, "arm_closure_summary-", run_time, ".csv")
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
    paste0(output_path, "arm_closures_sd", run_time, ".csv")
  )

  return(list(summ, sd_mx))
}


#' Plot predicted recruitment from file containing a CSV from
#' a single run
#' 
#' @param file_prefix Consistent beginning of filenames holding 
#' arm closure data. Defaults to `accrual`.
#' @param plot_prefix Prefix for file name to identify plot type. 
#' Defaults to `accrual_plot`.
#' @param run_time Specify a particular instance of `biomkrAccrual()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`. 
#' Used to select which files will be summarised.
#' @param output_path Directory where the output files from the 
#' `biomkrAccrual()` instance are located.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' 
#' @export
#' 
accrual_at_time_plot_from_file <- function(
  file_prefix = "accrual",
  plot_prefix = "all-accrual-from-file",
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

  makeifnot_dir(figs_path)

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

  # Plot 
  p <- plot(
    accrual_df, 
    plot_prefix = plot_prefix,
    run_time = run_time
  )

  print(p)

  ggplot2::ggsave(
    paste0(figs_path, plot_prefix, "-", run_time, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
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
#' @param x long format dataframe with columns "Week", 
#' "Arm" and "Recruitment".
#' @param ... Not used.
#' @param plot_prefix Prefix for file name to identify plot type.
#' Defaults to `accrual_plot`.
#' @param run_time Specify a particular instance of `biomkrAccrual()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`.
#' @param output_path = Directory where the output files from the 
#' `biomkrAccrual()` instance are located.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' @param target_arm_size Number of subjects required for each treatment arm.
#' @param target_control Number of subjects required for control arm(s).
#' @param target_interim Number of subjects required for treatment arm at 
#' interim analysis.
#' @param accrual_period Number of weeks in recruitment period.
#' @param interim_period Number of weeks to recruit for interim analysis.
#' 
#' @import ggplot2
#' @importFrom grDevices palette.colors
#' @importFrom rlang .data
#' 
#' @export
#' 
plot.accrualplotdata <- function(
  x,
  ...,
  plot_prefix = "accrual_plot",
  run_time = NULL,
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/"),
  target_arm_size = NA_integer_,
  target_control = NA_integer_,
  target_interim = NA_integer_,
  accrual_period = NA_integer_,
  interim_period = NA_integer_
) {
  
  accrual_df <- x
  arm_names <- levels(accrual_df$Arm)

  linetypes <- c(
    "Interim arm" = 2, "Experimental arm" = 3, "Control arm" = 4,
    "Interim accrual" = 5, "Total accrual" = 6
  )

  hline_y <- c(target_interim, target_arm_size, target_control)
  vline_x <- c(interim_period, accrual_period)

  p <- ggplot2::ggplot(
    accrual_df, 
    ggplot2::aes(
      x = .data$Week, 
      y = .data$Recruitment, 
      group = .data$Arm, 
      color = .data$Arm
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    # Use colourblind friendly Okabe-Ito palette
    ggplot2::scale_colour_manual(
      values = grDevices::palette.colors(length(arm_names))
    ) +
    ggplot2::geom_vline(
      xintercept = vline_x,
      linewidth = 1,
      linetype = 2:3,
      color = "grey75"
    ) +
    ggplot2::geom_hline(
      yintercept = hline_y,
      linewidth = 1,
      linetype = 4:6,
      color = "grey65"
    ) +
    ggplot2::labs(
      title = "Accrual plot"
    ) +
    theme_bma(base_size = 16)

  return(p)
}


#' Plot distributions of recruitment to arms at given time.
#' 
#' @param x Matrix with columns for each recruitment arm, 
#' including control.
#' @param ... For compliance with plot.default().
#' @param target Vector of targets for recruitment. First two
#' should be those directly relevant to the subject of the graph.
#' @param target_names Vector of target names, for labelling.
#' @param plot_id Type of plot for title, e.g. "Interim accrual".
#' @param adjust The adjust parameter from `ggplot2::geom_density`;
#' higher values mean more smoothing. Defaults to 1.
#' 
#' @importFrom stats reshape
#' @import ggplot2
#' @importFrom grDevices palette.colors
#' 
#' @export
#' 
plot.armtotals <- function(
  x,
  ...,
  target,
  target_names,
  plot_id,
  adjust = 1
) {
  data_df <- matrix_to_long(x)

  # Which of the accrual targets are within the dataset
  if (length(target) > 2) {
    target_indices <- 
      c(1:2, 2 + which(target[-c(1, 2)] <= max(data_df$Recruitment)))
  } else {
    target_indices <- seq_len(length(target))
  }
  target <- target[target_indices]
  target_names <- target_names[target_indices]

  p <- ggplot2::ggplot(
    data = data_df
  ) +
    ggplot2::geom_density(
      ggplot2::aes(
        x = Recruitment, group = Arm, fill = Arm, col = Arm
      ),
      alpha = 0.4, adjust = 1
    ) +
    ggplot2::scale_fill_manual(
      values = grDevices::palette.colors(length(unique(data_df$Arm))),
      aesthetics = c("color", "fill")
    ) +
    ggplot2::geom_vline(
      xintercept = target, 
      linetype = "dashed",
      linewidth = 1,
      colour = "grey65"
    ) +
    ggplot2::labs(
      y = "Probability density",
      title = plot_id,
    ) +
    theme_bma(base_size = 16)

  p <- label_vlines(p, target, target_names)

  return(p)
}


#' Adds labels for vlines to accrual plots
#' 
#' @param p Ggplot object.
#' @param target Vector of x axis positions of vlines.
#' @param target_names Vector of target names (excluding the
#' word `target`).
#' @param size Font size (in ggplot measure) for labels; 
#' defaults to 6.
#' 
label_vlines <- function(
  p,
  target,
  target_names,
  size = 6
) {
  # Get height of y axis for this particular plot
  label_y <- round(ggplot2::layer_scales(p)$y$range$range[2], 2)

  # Get x range of plot
  xrange <- round(ggplot2::layer_scales(p)$x$range$range, 2)

  whisker <- diff(xrange) * .2

  # There's more than one target.  Move labels for the one at the range end

  # Add labels for vlines
  abline_df <- data.frame(
    x = target, 
    y = label_y * 0.9, 
    label = paste(target_names, "\ntarget")
  )

  p <- p +
    ggplot2::geom_text(
      data = abline_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = size
    )

  return(p)
}


#' Plot single arm accrual plot
#' 
#' Bodge, fix later
#' 
#' @param data_df Dataframe of biomkrAccrual output.
#' @param arm_colours Vector of hexadecimal colours, one for each arm.
#' @param treatment_arms Vector of names of the treatment arms.
#' @param targets Vector of target names (excluding the
#' word `target`).
#' @param plot_id Vector of plot type names (typically "Interim" and "Accrual").
#' @param i Index of treatment arm to plot.
#' 
#' @import ggplot2
#' 
accrual_arm_plot <- function(
  data_df,
  arm_colours,
  treatment_arms,
  targets,
  plot_id,
  i
) {
  arm_names <- colnames(data_df)

  no_unique <- length(unique(data_df[, i]))

  if (no_unique == 1) {
    # BODGE - don't want to see this but need it to produce graph
    arm_col <- "white"
    alpha <- 0.0001
    unique_val <- unique(data_df[, i])
  } else {
    arm_col <- arm_colours[i]
    alpha <- 0.6
  }

  binwidth <- max(
    1,
    (max(data_df[, i]) - min(data_df[, i])) %/% 30
  )

  p <- ggplot2::ggplot(
    data = data_df
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(x = .data[[arm_names[i]]]),
      col = "white", fill = arm_col,
      alpha = alpha, 
      binwidth = binwidth
    ) 
  
  if (length(unique(data_df[, i])) == 1) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = unique(data_df[, i]),
        linewidth = 6,
        colour = arm_colours[i],
        alpha = 0.4
      ) 
  }

  p <- p + 
    ggplot2::geom_vline(
      xintercept = ifelse(
        treatment_arms[i], 
        targets[1], 
        targets[2]
      ), 
      linetype = "dashed",
      linewidth = 1,
      colour = "grey75"
    ) +
    ggplot2::labs(
      x = paste(
        "No. virtual patients recruited at", 
        "target week for",
        tolower(plot_id)
      ),
      y = "Probability density",
      title = paste(plot_id, "for", arm_names[i]),
    ) +
    theme_bma(base_size = 16)
  
  if (length(unique(data_df[, i])) == 1) {
    p <- p + ggplot2::scale_x_continuous(breaks = seq(
      unique_val - 1, length.out = 3
    ))
  } else {
    p <- p +
      ggplot2::scale_x_continuous(expand = expansion(mult = 0.07))
  }

  p <- label_vlines(
    p, 
    target = ifelse(
      treatment_arms[i],
      targets[1], 
      targets[2]
    ),
    target_names = ifelse(
      treatment_arms[i],
      plot_id, 
      paste0(plot_id, "\ncontrol")
    )
  )

  return(p)
}


#' Plot simulations of recruitment to given arm over time.
#' 
#' @param x Matrix with columns for each simulation.
#' @param ... For compliance with plot.default().
#' @param arm_colour Hexadecimal colour associated with arm.
#' @param target Vector of targets for recruitment. First two
#' should be those directly relevant to the subject of the graph.
#' @param target_names Vector of target names, for labelling.
#' @param plot_id Type of plot for title, e.g. "Treatment A".
#' 
#' @importFrom stats reshape
#' @import ggplot2
#' @importFrom grDevices palette.colors
#' 
#' @export
#' 
plot.armaccrual <- function(
  x,
  ...,
  arm_colour,
  target,
  target_names,
  plot_id
) {

  dimnames(x) <- list(
    Week = seq_len(nrow(x)),
    Simulation = c(seq_len(ncol(x)))
  )

  data_mx <- apply(x, 2, cumsum)


  ribbonwidth <- apply(data_mx, 1, function(x) quantile(x, c(0.025, 0.975)))


  data_df <- stats::reshape(
    as.data.frame(data_mx),
    direction = "long",
    varying = colnames(x),
    timevar = "Simulation",
    times = colnames(x),
    v.names = "Recruitment",
    idvar = "Week"
  )

  data_df$Simulation <- as.factor(data_df$Simulation)

  # Now summary data for ribbon
  ribbon_df <- as.data.frame(t(apply(
    data_mx, 1, function(x) quantile(x, c(0.025, 0.975))
  )))
  names(ribbon_df) <- c("Lower", "Upper")

  ribbon_df$Week <- seq_len(nrow(data_mx))
  ribbon_df$Mean <- rowMeans(data_mx)

  # Label hlines
  hline_df <- data.frame(
    x = nrow(data_mx) * .90,
    y = target,  
    label = paste(target_names, "\ntarget")
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = data_df,
      ggplot2::aes(
        x = Week,  
        y = Recruitment
      ),
      alpha = 0.2,
      col = "grey65"
    ) +
    ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(
        x = Week,
        ymin = Lower,
        ymax = Upper
      ),
      fill = arm_colour,
      alpha = 0.2
    ) +
    ggplot2::geom_hline(
      yintercept = target[1],
      linewidth = 1,
      linetype = 2,
      color = "grey75"
    ) +
    ggplot2::geom_hline(
      yintercept = target[2],
      linewidth = 1,
      linetype = 2,
      color = "grey75"
    ) +
    ggplot2::geom_text(
      data = hline_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = 6
    ) +
    ggplot2::geom_line(
      data = ribbon_df,
      ggplot2::aes(
        x = Week,
        y = Mean
      ),
      col = arm_colour,
      size = 1
    ) +
    ggtitle(paste(plot_id, "accrual")) +
    theme_bma(
      base_size = 14
    )
}



#' Extract from a simulation * week accrual matrix for a 
#' specific trial arm, the week in which an arm target 
#' threshold or thresholds are first met or exceeded for 
#' each simulation.
#' 
#' @param accrual Matrix of accrual data for a given arm
#' @param targets Target threshold(s) for recruitment.
#' 
#' @return List of vectors of the weeks in which each simulation
#' reaches the target accrual threshold.
#' 
threshold_week <- function(accrual, targets) {
  # Use cumulative version of accrual matrix
  accrual_cumsum_mx <- apply(
    accrual,
    2,
    cumsum
  )
  # Which week, if any, does simulation exceed target
  accrual_times_ls <- vector(mode = "list", length = length(targets))

  # Get list of vectors of weeks accrual meets threshold
  for (i in seq_len(length(targets))){
    accrual_times_ls[[i]] <- apply(
      accrual_cumsum_mx,
      2,
      function(x) {
        which(x >= targets[i])[1]
      }
    )
    # Set class so can use plot method
    accrual_times_ls[[i]] <- structure(
      accrual_times_ls[[i]],
      class = c("targetweek", "integer")
    )

  }
  return(accrual_times_ls)
}


#' Plot predicted week at which recruitment reaches specified targets
#' for each arm from file containing a json consisting of a list
#' of matrices of arm recruitment by week, one for each arm.
#' 
#' @param file_prefix Consistent beginning of filenames holding 
#' arm closure data. Defaults to `accrual`.
#' @param plot_prefix Prefix for file name to identify plot type. 
#' Defaults to `accrual_plot`.
#' @param run_time Specify a particular instance of `biomkrAccrual()`
#' execution using a date-time format `yyyy-mm-dd-hh-mm-ss`. 
#' Used to select which files will be summarised.
#' @param output_path Directory where the output files from the 
#' `biomkrAccrual()` instance are located.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' @param target_treatment Recruitment target for treatment arms.
#' @param target_control Recruitment target for control arms. 
#' @param arm_colours Vector of hexadecimal colours, one for each arm.
#' If only one value is supplied, it will be used for all arms. If no
#' value is supplied, a colourblind-friendly palette will be used.
#' @export
#' 
accrual_to_target_plot_from_file <- function(
  file_prefix = "arm_accrual",
  plot_prefix = "accrual-arm-target-week-from-file",
  run_time = "26-02-08_14-15-05",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/"),
  target_treatment = NA_integer_,
  target_control = NA_integer,
  arm_colours = NULL
) {
  # Validate input

  checkmate::assert_directory_exists(
    file.path(output_path), 
    access = "rx"
  )

  input_file <- paste0(
    output_path, file_prefix, "_", run_time, ".json"
  )

  checkmate::assert_file_exists(file.path(input_file))

  # Recruitment targets should be integerish and not missing
  checkmate::assert_integerish(
    target_treatment, lower = 1, upper = 10^7, len = 1, any.missing = FALSE
  )
  checkmate::assert_integerish(
    target_control, lower = 1, upper = 10^7, len = 1, any.missing = FALSE
  )

  # Colours should be hexadecimal
  if (!is.null(arm_colours)) {
    if (!all(grepl("^#[A-Fa-f0-9]+$", arm_colours))) {
      rlang::abort("Colours should be in the form #123456.")
    }
  }

  makeifnot_dir(figs_path)

  accrual_raw_ls <- jsonlite::read_json(
    input_file,
    simplifyVector = TRUE
  )

  plot_id <- names(accrual_raw_ls)
  target_index <- 2 - startsWith(plot_id, "T")
  targets <- c(target_treatment, target_control)

  if (is.null(arm_colours)) {
    col_order <- c(seq_len(length(plot_id))[-1], 1)
    palette <- grDevices::palette.colors(
      palette = "R4",
      length(plot_id) + 1
      # One pink is more than enough
    )[-6]
    arm_colours <- palette[col_order]
  } else if (length(arm_colours) == 1) {
    arm_colours <- rep(arm_colours, length(plot_id))
  } else if (length(arm_colours) != length(plot_id)) {
    rlang::abort(paste(
      "Invalid number of arm_colours: please supply",
      "either nothing, one value or one value per arm."
    ))
  }

  for (i in seq_len(length(plot_id))) {
    accrual_times <- threshold_week(
      accrual = accrual_raw_ls[[i]],
      targets = targets[target_index[i]]
    )
    
    p <- plot(
      accrual_times[[1]],
      arm_colour = arm_colours[i],
      target = targets[target_index[i]],
      plot_id = plot_id[i]
    )

    print(p)
    
  }



}


#' Plot distribution of time to accrual to a specified
#' target for a specified arm.
#' 
#' @param x Vector of times to accrual.
#' @param ... For compliance with plot.default().
#' @param arm_colour Hexadecimal colour associated with arm.
#' @param target Accrual target, for labelling.
#' @param target_names Target name, for labelling.
#' @param plot_id Type of plot for title, e.g. "Treatment A".
#' 
#' @importFrom stats reshape
#' @import ggplot2
#' @importFrom grDevices palette.colors
#' 
#' @export
#' 
plot.targetweek <- function(
  x,
  ...,
  arm_colour,
  target,
  target_names = NULL,
  plot_id
) {
  if (all(is.na(x))) {
    message(paste(
      "No simulations of arm", plot_id, 
      "reached the target", target
    ))
    return(NULL)
  }

  target_names <- tolower(target_names)

  plot_label <- paste0(
    "Week in which accrual to ",
    plot_id,
    " met the ",
    ifelse(is.null(target_names), "", paste0(target_names, " ")),
    "target of ", 
    target
  )
  plot_sublabel <- paste(
    sum(is.na(x)),
    "of",
    length(x),
    "simulations did not reach the target before all sites were closed"
  )

  no_unique <- length(unique(x[!is.na(x)]))

  if (no_unique == 1) {
    # BODGE - don't want to see this but need it to produce graph
    arm_col <- "white"
    alpha <- 0.0001
    unique_val <- unique(x[!is.na(x)])
  } else {
    arm_col <- arm_colour
    alpha <- 0.6
  }


  # Pretty breaks for x axis
  xrange <- range(x, na.rm = TRUE)
  ticksep <- 1 + diff(xrange) %/% 20
  xmin <- xrange[1] - (xrange[1] %% ticksep)
  xmax <- xrange[2] + (ifelse(xrange[2] %% ticksep > 0, ticksep, 0))

  # Pretty breaks for y axis

  p <- ggplot2::ggplot(data.frame(x = x[!is.na(x)])) +
    ggplot2::geom_histogram(
      ggplot2::aes(x = x),
      col = "white", fill = arm_col,
      alpha = alpha,
      binwidth = 1
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(
        xmin, 
        xmax, 
        ticksep
      )
    )
  
  if (no_unique == 1) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = unique_val,
        linewidth = 6,
        colour = arm_colour,
        alpha = 0.4
      )
  }
    
    
  p <- p +
    labs(
      title = plot_label,
      subtitle = plot_sublabel,
      x = "Week",
      y = "Count"
    ) +
    theme_bma(base_size = 14)

  return(p)

}