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
#' @param accrual_period Maximum number of weeks to recruit 
#' (defaults to 80)
#' @param target_arm_size Maximum size for all experimental arms
#' (defaults to 308)
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
  output_path = "biomkrAccrual_output_data/",
  accrual_period = 10,
  target_arm_size = 100,
  shared_control = TRUE,
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE,
  fixed_region_prevalences = FALSE
) {

  # Information for setting up dataframe
  arms_ls <- 
    jsonlite::read_json(paste0(data_path, arms_file), simplifyVector = TRUE)

  # Define matrix of zeroes for efficiency
  arm_closures_mx <- matrix(0, nrow = n, ncol = length(arms_ls))
  arm_totals_mx <- matrix(0, nrow = n, ncol = ifelse(
    shared_control,
    length(arms_ls) + 1,
    2 * length(arms_ls)
  ))

  # Set column names
  colnames(arm_closures_mx) <- names(arms_ls)
  if (shared_control) {
    colnames(arm_totals_mx) <- c(names(arms_ls), "Control")
  } else {
    colnames(arm_totals_mx) <- c(names(arms_ls), paste0("C-", names(arms_ls)))
  }

  # Run batches
  for (irun in seq(n)) {
    out_ls <- spine(
      centres_file = centres_file,
      arms_file = arms_file,
      data_path = data_path,
      accrual_period = accrual_period,
      target_arm_size = target_arm_size,
      shared_control = shared_control,
      fixed_centre_starts = fixed_centre_starts,
      fixed_site_rates = fixed_site_rates,
      fixed_region_prevalences = fixed_region_prevalences
    )

    arm_closures_mx[irun, ] <- out_ls[[1]]
    arm_totals_mx[irun, ] <- out_ls[[2]]
  }

  # Keep copies of output, stamped with datetime
  datetime <- format(Sys.time(), "%y-%m-%d_%H-%M-%S")

  write.csv(
    as.data.frame(arm_closures_mx),
    paste0(output_path, "arm_closures_", datetime, ".csv")
  )

  write.csv(
    as.data.frame(arm_totals_mx),
    paste0(output_path, "arm_totals_", datetime, ".csv")
  )

  print(summary(as.data.frame(arm_closures_mx)))
  print(summary(as.data.frame(arm_totals_mx)))
}