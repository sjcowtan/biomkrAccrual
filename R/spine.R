#' @include prevalence.R

#' Driver for the procedure
#' @param target_arm_size Number of patients required per treatment arm
#' 
#' @examples 
#' spine()
#' 
#' @export
spine <- function(
  target_arm_size = 308,
  target_interim = target_arm_size / 2,
  target_control = 704,
  margin = 3, 
  # Specify this for equal rates
  av_site_rate_month = 2,
  accrual_period = 36,
  # active : control ratio (all active the same)
  ctrl_ratio = c(1, 1),
  no_samples = 100,
  centres_file = "centres.csv",
  prop_file = "proportions.csv",
  centre_start_file = "centre_starts.csv",
  arms_file = "arms.json",
  phase_file = "phase_change_weeks.csv",
  average_file = "mean_recruitment.csv",
  # Use this if expected site rates not equal
  site_rate_file = "site_rates.csv",
  data_path = "inst/extdata/",
  figs_path = "figures/",
  fixed_centre_starts = TRUE,
  fixed_site_rates = TRUE
  
) {
  # Verify inputs
  ## fail if directory doesn't exist
  ## append "/" if no slash on end
  ## If _file doesn't end in .csv add it 
  ## Check if files exist and are readable
  ## Check for switches e.g. av_site_rate_month first

  # Read parameters
  prop_params_df <- read.csv(paste0(data_path, prop_file))
  arms_ls <- jsonlite::read_json(paste0(data_path, arms_file))
  goals_df <- read.csv(paste0(data_path, centres_file))

  # Fail if read fails
  # Fail if wrong format
  if (goals_df[nrow(goals_df), ]$no_centres != 0) {
    rlang::abort(paste0(
      "The last line of ", data_path, centres_file,
      " should consist of the pair (month after site opening finishes, 0)"
    ))
  }

  # Create structure object
  trial_structure_instance <- trial_structure(prop_params_df, arms_ls)
 
  print("Before do_sample")
  print(trial_structure_instance@treatment_arm_struct)
  print(trial_structure_instance@treatment_arm_prev)

  treats_out_ls <- do_sample(
      prop_params_df,
      arms_ls,
      trial_structure_instance,
      target_arm_size,
      margin,
      av_site_rate_month,
      accrual_period,
      goals_df,
      ctrl_ratio,
      data_path = data_path,
      centre_start_file = centre_start_file,
      fixed_centre_starts = fixed_centre_starts,
      fixed_site_rates = fixed_site_rates
    )
    print("After do_sample")
    print(trial_structure_instance@treatment_arm_struct)
    print(trial_structure_instance@treatment_arm_prev)

}


#' Driver for single instance of sampling
#' 
do_sample <- function(
  prop_params_df, 
  arms_ls, 
  trial_structure_instance, 
  target_arm_size, 
  margin, 
  av_site_rate_month, 
  accrual_period, 
  goals_df, 
  ctrl_ratio, 
  data_path, 
  centre_start_file, 
  fixed_centre_starts, 
  fixed_site_rates
) {
  

  arms_to_remove <- as.integer(c(1, 4))
  trial_structure_instance <- 
    remove_treat_arms(trial_structure_instance, arms = arms_to_remove)


  return(trial_structure_instance)

}