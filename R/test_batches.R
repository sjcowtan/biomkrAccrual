#' Run a number of batches of recruitment prediction and
#' collect summary statistics on arm closures, and final
#' recruitment totals for all experimental and control arms
#' @param n Number of instances to run (defaults to 100)
#' @param centres_file File containing recruitment centre data
#' (defaults to "centres.csv")
#' @param arms_file File containing list of which recruitment
#' arms recruit to which experimental arm
#' @param accrual_period Maximum number of weeks to recruit 
#' (defaults to 80)
#' @param target_arm_size Maximum size for all experimental arms
#' (defaults to 308)
#' @param shared_control = TRUE,
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
  data_path = "inst/extdata/",
  output_path = "output_data/",
  accrual_period = 10,
  target_arm_size = 100,
  shared_control = TRUE,
  fixed_site_rates = FALSE
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
      shared_control = shared_control
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