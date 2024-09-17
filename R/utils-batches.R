#' Run a number of batches of recruitment prediction and
#' collect summary statistics on arm closures, and final
#' recruitment totals for all experimental and control arms
#' @param n Number of instances to run (defaults to 100)
#' @param centres_file File containing recruitment centre data
#' (defaults to "centres.csv")
#' @param arms_file File containing list of which recruitment
#' arms recruit to which experimental arm
#' @param data_path Folder where `centres_file`, `prop_file` and
#' `arms_file` are located. Defaults to the location of the package
#' example data in the package installation; this should be changed. 
#' @param output_path Folder where data generated during execution
#' will be stored; defaults to `../biomkrAccrual_output_data/`.
#' @param figs_path Folder where figures generated during execution
#' will be stored; defaults to the `figures` subdirectory in
#' `output_path`.
#' @param accrual_period Maximum number of weeks to recruit 
#' (defaults to 80)
#' @param target_arm_size Maximum size for all experimental arms
#' (defaults to 308)
#' @param target_interim Recruitment target for experimental arms 
#' at interim analysis; defaults to `target_arm_size / 2`.
#' @param target_control Maximum size for all control arms
#' (defaults to 308)
#' @param target_interim_control Recruitment target for control arms 
#' at interim analysis; defaults to `target_control / 2`.
#' @param shared_control = TRUE,
#' @param fixed_centre_starts TRUE if centres are assumed to start
#' exactly when planned; FALSE if some randomisation should be added.
#' @param fixed_site_rates TRUE if centre recruitment rates should 
#' be treated as exact; FALSE if they should be drawn from a gamma
#' distribution with a mean of the specified rate.
#' @param fixed_region_prevalences TRUE if biomarker prevalences 
#' should be considered to be identical for all sites within a 
#' region; FALSE if they should be drawn from a Dirichlet distribution
#' with a mean of the specified prevalence.
#' 
#' @return Dataframe of site closing times
#' @return Dataframe of experimental arm totals
#' @export
#' 
#' @importFrom jsonlite read_json
#' @importFrom utils write.csv
#' 
batch <- function(
  n = 100,
  centres_file = "centres.csv",
  arms_file = "arms.json",
  data_path = "extdata/",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/"),
  accrual_period = 6,
  interim_period = 3,
  target_arm_size = 60,
  target_interim = target_arm_size / 2,
  target_control = 180,
  target_interim_control = target_control / 2,
  shared_control = TRUE,
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE,
  fixed_region_prevalences = FALSE
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
      centres_file = centres_file,
      arms_file = arms_file,
      data_path = data_path,
      accrual_period = accrual_period,
      interim_period = interim_period,
      target_arm_size = target_arm_size,
      target_interim = target_interim,
      target_control = target_control,
      shared_control = shared_control,
      fixed_centre_starts = fixed_centre_starts,
      fixed_site_rates = fixed_site_rates,
      fixed_region_prevalences = fixed_region_prevalences,
      quietly = TRUE
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
        seq_len(accrual_instance@interim_period), , 
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

  p <- plot(
    arm_interim_mx, 
    target = c(
      target_interim, target_interim_control, 
      target_arm_size, target_control
    ), 
    target_names = c(
      "Interim", "Interim control", 
      "Accrual", "Accrual control"
    )
  )

  ggplot2::ggsave(
    paste0(figs_path, "arm-totals-interim-", run_time, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
  )

  print(p)
}