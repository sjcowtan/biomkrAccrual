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
  shared_control = TRUE,
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
  arms_ls <- 
    jsonlite::read_json(paste0(data_path, arms_file), simplifyVector = TRUE)
  centres_df <- read.csv(paste0(data_path, centres_file))

  # Add fail if read fails

  # Fail if centres_file is in wrong format
  if (isFALSE(all.equal(
    names(centres_df), 
    c("site", "start_month", "mean_rate", "prevalence_set", "site_cap")
    # site_cap is optional, no cap if not present
  )) || isFALSE(all.equal(
    names(centres_df), c("site", "start_month", "mean_rate", "prevalence_set")
  ))) {
    rlang::abort(paste(
      "Format error: centres.csv should have columns site,",
      "start_month, mean_rate, prevalence_set, and optionally site_cap"
    ))
  }

  # Get start weeks & order centres_df by start week and site number
  centres_df <- do_clean_centres(centres_df)

  print(centres_df)

  # Make control ratio 1:x if used
  if (!is.null(ctrl_ratio) && !identical(ctrl_ratio[1], 1)) {
    ctrl_ratio <- ctrl_ratio / ctrl_ratio[1]
  } else if (is.null(ctrl_ratio) && !is.null(target_control)) {
    ctrl_ratio <- c(1, target_control / target_arm_size)
  }

  # Generate target_control if needed
  if (is.null(target_control) && shared_control) {
    if (is.null(ctrl_ratio)) {
      rlang::abort(paste(
        "For shared control, either ctrl_ratio or", 
        "target_control must be specified."
      ))
    } else {
      target_control <- target_arm_size * ctrl_ratio[2]
    }
  } 

  # Total target recruitment
  target_recruit <- ifelse(
    shared_control, 
    target_arm_size * length(arms_ls) + target_control,
    target_arm_size * length(arms_ls) * (2 * ctrl_ratio[2])
  )

  # Complete site cap if incomplete, using recruitment target
  if (!("site_cap" %in% names(centres_df))) {
    centres_df$site_cap <- target_recruit
  } else {
    centres_df$site_cap[is.na(centres_df$site_cap)] <- target_recruit
  }

  # Create structure object
  trial_structure_instance <- 
    trial_structure(prop_params_df, arms_ls, centres_df, shared_control)
  print("Trial structure")
  print(trial_structure_instance@treatment_arm_struct)
  print("Experimental arm prev")
  print(trial_structure_instance@experimental_arm_prevalence)


  print(c("Accrual period", accrual_period))
  print(
    length(trial_structure_instance@treatment_arm_ids) + 
      ifelse(
        shared_control, 
        1, 
        length(trial_structure_instance@treatment_arm_ids)
      )
  )

  # Create accrual object
  accrual_instance <- accrual(
    treatment_arm_ids = trial_structure_instance@treatment_arm_ids,
    shared_control = shared_control,
    centres_df = centres_df,
    accrual_period = accrual_period
  )

  
  accrual_instance <- accrue_week(
    accrual_instance, 
    trial_structure_instance,
    target_arm_size
  )
  print("Accrual week 1")
  print(head(accrual_instance@accrual))

  return()

  print(paste("Next week", accrual_instance@week))

  accrual_instance <- accrue_week(accrual_instance, target_arm_size)

  print("Accrual week 2")
  print(head(accrual_instance@accrual))
  print(paste("Next week", accrual_instance@week))

  accrual_instance <- accrue_week(
    list(accrual_instance, trial_structure_instance), target_arm_size
  )

  print("Accrual week 3")
  print(head(accrual_instance@accrual))
  print(paste("Next week", accrual_instance@week))

  accrual_instance <- accrue_week(accrual_instance, target_arm_size)
  
  print(c("Site sums", site_sums(accrual_instance)))
  print(c(
    "Arm sums", treat_sums(accrual_instance)[
      seq_len(length(accrual_instance@phase_changes))
    ]
  ))
}
