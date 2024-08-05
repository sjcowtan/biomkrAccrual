#' @include prevalence.R



#' @title Driver for the procedure
#' @param target_arm_size Number of patients required per treatment arm
#' 
#' @examples 
#' spine()
#' 
#' @importFrom jsonlite read_json
#' @importFrom rlang abort
#' @importFrom utils read.csv
#' 
#' @export
spine <- function(
  target_arm_size = 308,
  target_interim = target_arm_size / 2,
  target_control = 704,
  accrual_period = 36,
  precision = 10,
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
  data_path = "extdata/",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "figures/"),
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE
) {
  # Verify inputs
  ## append "/" if no slash on end
  #data_path <- gsub("(\\w+)$", "\\1/", data_path)
  output_path <- gsub("(\\w+)$", "\\1/", output_path)
  # Make into full path so only one set of syntax needed
  if (!grepl("^/", output_path)) {
    output_path <- paste0(getwd(), "/", output_path)
  }
  
  ## If _file doesn't end in .csv add it 
  

  ## Check for switches e.g. av_site_rate_month first
  checkmate::assert_directory_exists(system.file(
    data_path, package = "biomkrAccrual"
  ), access = "rx")
  checkmate::assert_file_exists(system.file(
    data_path, prop_file, package = "biomkrAccrual"
  ), access = "r")
  checkmate::assert_file_exists(system.file(
    data_path, arms_file, package = "biomkrAccrual"
  ), access = "r")
  checkmate::assert_file_exists(system.file(
    data_path, centres_file, package = "biomkrAccrual"
  ), access = "r")

  # Set up output directory if does not already exist
  if (checkmate::test_directory_exists(file.path(
    output_path
  ))) {
    # Can we write to it?
    checkmate::assert_access(file.path(
      output_path, access = "wx"
    )) 
  } else {
    # Wrap me in a tryCatch
    tryCatch(
      dir.create(output_path),
      # Convert warning to error
      warning = function(w) rlang::abort(w$message)
    )
  }

  # Set up output figures directory if does not already exist
  if (checkmate::test_directory_exists(file.path(
    figs_path
  ))) {
    # Can we write to it?
    checkmate::assert_access(file.path(
      figs_path, access = "wx"
    )) 
  } else {
    # Wrap me in a tryCatch
    tryCatch(
      dir.create(figs_path),
      # Convert warning to error
      warning = function(w) rlang::abort(w$message)
    )
  }
  
  # Read parameters
  prop_params_df <- read.csv(system.file(
    data_path, prop_file, package = "biomkrAccrual"
  )) 
    
  arms_ls <- 
    jsonlite::read_json(system.file(
      data_path, arms_file, package = "biomkrAccrual"
    ), simplifyVector = TRUE)
  
  centres_df <- read.csv(system.file(
    data_path, centres_file, package = "biomkrAccrual"
  ))

  # Add fail if read fails

  # Fail if centres_file is in wrong format
  if (isFALSE(all.equal(
    names(centres_df), 
    c("site", "start_month", "mean_rate", "region", "site_cap")
    # site_cap is optional, no cap if not present
  )) || isFALSE(all.equal(
    names(centres_df), c("site", "start_month", "mean_rate", "region")
  ))) {
    rlang::abort(paste(
      "Format error: centres.csv should have columns site,",
      "start_month, mean_rate, region, and optionally site_cap"
    ))
  }

  # Get start weeks & order centres_df by start week and site number
  centres_df <- do_clean_centres(centres_df)
  centres_df$start_week <- get_weeks(centres_df$start_month - 1) + 1

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
    trial_structure(
      prop_params_df, 
      arms_ls, 
      centres_df, 
      precision, 
      shared_control
    )

  # Create accrual object
  accrual_instance <- accrual(
    treatment_arm_ids = trial_structure_instance@treatment_arm_ids,
    shared_control = shared_control,
    centres_df = centres_df,
    accrual_period = get_weeks(accrual_period)
  )

  while (
    # Any arms are recruiting
    any(trial_structure_instance@treatment_arm_struct) &&
      # Any sites are recruiting
      length(accrual_instance@active_sites) > 0 &&
      # Not out of time
      accrual_instance@week <= accrual_instance@accrual_period
  ) {

    # Add a week's accrual
    obj_list <- accrue_week(
      accrual_instance, 
      trial_structure_instance,
      target_arm_size,
      fixed_site_rates
    )

    accrual_instance <- obj_list[[1]]
    trial_structure_instance <- obj_list[[2]]

    # Increment pointer for the next week to accrue
    accrual_instance@week <- accrual_instance@week + as.integer(1)
  }

  # Trim accrual to actual recruitment length
  accrual_instance@accrual <- 
    accrual_instance@accrual[seq(accrual_instance@week - 1), , ]

  # Return summary statistics
  return(list(
    accrual_instance@phase_changes, 
    treat_sums(accrual_instance)
  ))
}
