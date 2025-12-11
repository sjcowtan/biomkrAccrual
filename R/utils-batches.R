#' Run a number of batches of recruitment prediction and
#' collect summary statistics on arm closures, and final
#' recruitment totals for all experimental and control arms
#' @param n Number of instances to run (defaults to 100)
#' @param target_arm_size Number of patients required per 
#' treatment arm
#' @param target_interim Number of patients required per 
#' arm for interim analysis; defaults to `target_arm_size \ 2`
#' @param target_control Number of patients required for the
#' control arm(s)
#' @param shared_control TRUE if all treatment arms share the
#' same control arm; FALSE if each treatment arm has its own 
#' control. Defaults to TRUE.
#' @param accrual_period Recruitment period (months).
#' @param interim_period Recruitment period to interim (months).
#' @param precision For the Dirichlet model of biomarker prevalences, 
#' variability decreases as precision increases. Defaults to 10.
#' @param var_lambda Variance estimate for site recruitment rates.  
#' Defaults to 0.25.
#' @param control_ratio Ratio of patient allocation to treatment arm
#' versus control for all active arms; defaults to c(1, 1).
#' @param centres_file Name of CSV file with information about 
#' each recruitment centre; this should have columns "site", 
#' "start_month", "mean_rate", "region" and optionally "site_cap"
#' if recruitment is capped per site. Defaults to `centres.csv`.
#' @param prop_file Name of CSV file with expected biomarker prevalences
#' for each region; this should have one column "category", naming
#' the biomarkers or biomarker combinations, and one column per
#' region. Defaults to `proportions.csv`.
#' @param arms_file Name of JSON file describing which recruitment
#' arms (defined by biomarkers) recruit to which treatment arms. 
#' Defaults to `arms_json`.
#' @param data_path Folder where `centres_file`, `prop_file` and
#' `arms_file` are located. Defaults to the location of the package
#' example data in the package installation; this should be changed. 
#' @param output_path Folder where data generated during execution
#' will be stored; defaults to `../biomkrAccrual_output_data/`.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' @param fixed_centre_starts TRUE if centres are assumed to start
#' exactly when planned; FALSE if some randomisation should be added.
#' @param fixed_site_rates TRUE if centre recruitment rates should 
#' be treated as exact; FALSE if they should be drawn from a gamma
#' distribution with a mean of the specified rate.
#' @param fixed_region_prevalences TRUE if biomarker prevalences 
#' should be considered to be identical for all sites within a 
#' region; FALSE if they should be drawn from a Dirichlet distribution
#' with a mean of the specified prevalence.
#' @param quietly Defaults to FALSE, which displays the output from
#' each run. Set to TRUE to generate data and figures without displaying
#' them.
#' @param keep_files Save data files and plots generated during the run. 
#' Defaults to TRUE.
#' @return Dataframe of site closing times
#' @return Dataframe of experimental arm totals
#' @export
#' 
#' @importFrom jsonlite read_json
#' @importFrom utils write.csv
#' 
biomkrAccrualSim <- function(
  n = 100,
  target_arm_size = 60,
  target_interim = target_arm_size / 2,
  target_control = 180,
  shared_control = TRUE,
  accrual_period = 50 / 4,
  interim_period = accrual_period / 2,
  precision = 10,
  var_lambda = 0.25,
  # active : control ratio (all active the same)
  control_ratio = c(1, 1),
  centres_file = "centres.csv",
  prop_file = "proportions.csv",
  arms_file = "arms.json",
  data_path = "extdata/",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/"),
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE,
  fixed_region_prevalences = FALSE,
  quietly = FALSE,
  keep_files = TRUE
) {
  # Timestamp for batch files (but not individual run files)
  run_time <- format(Sys.time(), "%F-%H-%M-%S")

  # Information for setting up dataframe
  arms_ls <- 
    jsonlite::read_json(system.file(
      data_path, arms_file, package = "biomkrAccrual"
    ), simplifyVector = TRUE)

  # Define matrix of zeroes for efficiency
  arm_closures_mx <- structure(
    matrix(0, nrow = n, ncol = length(arms_ls)),
    class = c("armtotals", "matrix", "array")
  )
  arm_totals_mx <- structure(
    matrix(0, nrow = n, ncol = ifelse(
      shared_control,
      length(arms_ls) + 1,
      2 * length(arms_ls)
    )),
    class = c("armtotals", "matrix", "array")
  )
  arm_interim_mx <- structure(
    matrix(0, nrow = n, ncol = ifelse(
      shared_control,
      length(arms_ls) + 1,
      2 * length(arms_ls)
    )),
    class = c("armtotals", "matrix", "array")
  )

  # Set column names
  colnames(arm_closures_mx) <- names(arms_ls)
  if (shared_control) {
    colnames(arm_totals_mx) <- c(names(arms_ls), "Control")
  } else {
    colnames(arm_totals_mx) <- c(names(arms_ls), paste0("C-", names(arms_ls)))
  }
  colnames(arm_interim_mx) <- colnames(arm_totals_mx)

  # Run batches
  for (irun in seq(n)) {
    accrual_instance <- biomkrAccrual(
      target_arm_size = target_arm_size,
      target_interim = target_interim,
      target_control = target_control,
      shared_control = shared_control,
      accrual_period = accrual_period,
      interim_period = interim_period,
      precision = precision,
      var_lambda = var_lambda,
      control_ratio = control_ratio,
      centres_file = centres_file,
      arms_file = arms_file,
      data_path = data_path,
      fixed_centre_starts = fixed_centre_starts,
      fixed_site_rates = fixed_site_rates,
      fixed_region_prevalences = fixed_region_prevalences,
      quietly = TRUE,
      keep_files = FALSE
    )

    arm_closures_mx[irun, ] <- accrual_instance@phase_changes
    arm_totals_mx[irun, ] <- treat_sums(
      accrual_instance@accrual[
        seq_len(min(
          nrow(accrual_instance@accrual), accrual_instance@accrual_period
        )), ,
      ]

    )
    arm_interim_mx[irun, ] <- treat_sums(
      accrual_instance@accrual[
        seq(accrual_instance@interim_period), , 
      ]
    )
  }

  # Keep copies of output, stamped with datetime
  datetime <- format(Sys.time(), "%y-%m-%d_%H-%M-%S")

  write.csv(
    as.data.frame(arm_closures_mx),
    paste0(output_path, "arm_closures_", datetime, ".csv"),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(arm_totals_mx),
    paste0(output_path, "arm_totals_", datetime, ".csv"),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(arm_interim_mx),
    paste0(output_path, "arm_interim_totals_", datetime, ".csv"),
    row.names = FALSE
  )

  print(summary(as.data.frame(arm_closures_mx)))
  print(summary(as.data.frame(arm_totals_mx)))
  print(summary(as.data.frame(arm_interim_mx)))

  # Interim plot

  ### This is a bodge
  target_interim_control <- target_control * target_interim / target_arm_size

  p <- plot(
    arm_interim_mx, 
    target = c(
      target_interim, target_interim_control, 
      target_arm_size, target_control
    ), 
    target_names = c(
      "Interim", "Interim\ncontrol", 
      "Accrual", "Accrual\ncontrol"
    ),
    target_week = interim_period
  )

  ggplot2::ggsave(
    paste0(figs_path, "arm-totals-interim-", run_time, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
  )

  print(p)

  #


  # Individual accrual plots
  arm_names <- dimnames(arm_totals_mx)[[2]]

  ## Mark arms as treatment or control
  treatment_arms <- startsWith(arm_names, "T")

  ## Same colours as in interim plot
  col_order <- c(seq_len(length(treatment_arms))[-1], 1)
  arm_colours <- grDevices::palette.colors(length(treatment_arms))[col_order]
  
  

  # Total accrual plots

  data_ls <- list(
    Interim = as.data.frame(arm_interim_mx),
    Accrual = as.data.frame(arm_totals_mx)
  )
  target_ls <- list(
    Interim = c(target_interim, target_interim_control),
    Accrual = c(target_arm_size, target_control)
  )

  ## Loop across interim and total
  for (j in seq_len(length(data_ls))) {
    # Loop across all arms
    for (i in seq(treatment_arms)) {
      p <- accrual_arm_plot(
        data_ls[[j]],
        arm_colours, 
        treatment_arms,
        target_ls[[j]],
        plot_id = names(data_ls)[j],
        i
      )

      ggplot2::ggsave(
        paste0(
          figs_path, 
          "arm-totals-",
          tolower(names(data_ls)[j]),
          "-",
          arm_names[i], "-", 
          run_time, 
          ".png"
        ),
        plot = p,
        width = 12,
        height = 8,
        dpi = 400
      )

      print(p)
    }
  }
}


#' Format batch accrual data in long format for plotting.
#' 
#' @param data Matrix of accrual data.
#' 
matrix_to_long <- function(data) {

  arm_names <- dimnames(data)[[2]]

  data_df <- stats::reshape(
    as.data.frame(data),
    direction = "long",
    varying = arm_names,
    timevar = "Arm",
    times = arm_names,
    v.names = "Recruitment",
    idvar = "Run"
  )

  return(data_df)
}