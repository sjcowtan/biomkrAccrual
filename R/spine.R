#' @include prevalence.R



#' @title Command line 
#' 
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
#' @param ctrl_ratio Ratio of patient allocation to treatment arm
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
#' 
#' @examples 
#' biomkrAccrual()
#' 
#' @import checkmate 
#' @importFrom jsonlite read_json
#' @importFrom rlang abort warn
#' @importFrom utils read.csv
#' @importFrom ggplot2 ggsave
#' 
#' @export
biomkrAccrual <- function(
  target_arm_size = 60,
  target_interim = target_arm_size / 2,
  target_control = 180,
  shared_control = TRUE,
  accrual_period = 50 / 4,
  interim_period = accrual_period / 2,
  precision = 10,
  # active : control ratio (all active the same)
  ctrl_ratio = c(1, 1),
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

  checkmate::assert_logical(
    fixed_region_prevalences,
    any.missing = FALSE,
    null.ok = FALSE
  )

  if (fixed_region_prevalences && 
    checkmate::test_numeric(
      precision, 
      any.missing = FALSE, 
      null.ok = FALSE
    )
  )  {
    rlang::warn(paste("Value given for precision when", 
      "fixed_region_prevalences is TRUE: fixed_region_prevalences",
      "will take precendence."
    ))
    precision <- NULL
  }

  checkmate::assert_numeric(
    precision,
    any.missing = FALSE,
    lower = 10^-6,
    finite = TRUE,
    len = 1,
    null.ok = TRUE
  )

  if (!fixed_region_prevalences && is.null(precision)) {
    rlang::abort(paste("Either fixed_region_prevalences", 
      "must be TRUE, or a value must be given for the",
      "precision of the Dirichlet distribution."
    ))
  }


  # Verify inputs
  ## append "/" if no slash on end
  data_path <- gsub("(\\w+)$", "\\1/", data_path)
  output_path <- gsub("(\\w+)$", "\\1/", output_path)
  figs_path <- gsub("(\\w+)$", "\\1/", figs_path)


  # Make into full path so only one set of syntax needed
  if (!grepl("^/", output_path)) {
    output_path <- paste0(getwd(), "/", output_path)
  }
  
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
  makeifnot_dir(output_path, min_access = "rwx")
  
  # Set up output figures directory if does not already exist
  makeifnot_dir(figs_path, min_access = "rwx")
  
  # Read parameters
  prop_params_df <- utils::read.csv(system.file(
    data_path, prop_file, package = "biomkrAccrual"
  )) 
    
  arms_ls <- 
    jsonlite::read_json(system.file(
      data_path, arms_file, package = "biomkrAccrual"
    ), simplifyVector = TRUE)
  
  centres_df <- utils::read.csv(system.file(
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

  # Set run_time to timestamp output files
  run_time <- format(Sys.time(), "%F-%H-%M-%S")

  # Get start weeks & order centres_df by start week and site number
  centres_df <- do_clean_centres(centres_df)
  centres_df$start_week <- get_weeks(centres_df$start_month - 1) + 1

  # Make control ratio sum to 1
  if (is.null(ctrl_ratio)) {
    if (!is.null(target_control)) {
      ctrl_ratio <- c(1, target_control / target_arm_size)
    } else {
      rlang::abort(paste(
        "For shared control, either ctrl_ratio or", 
        "target_control must be specified."
      ))
    }
  }
  ctrl_ratio <- ctrl_ratio / sum(ctrl_ratio)

  # Generate target_control if needed
  if (is.null(target_control) && shared_control) {
    target_control <- target_arm_size * ctrl_ratio[2]
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
      props_df = prop_params_df, 
      arms_ls = arms_ls, 
      centres_df = centres_df, 
      precision = precision, 
      shared_control = shared_control,
      ctrl_ratio = ctrl_ratio,
      fixed_region_prevalences = fixed_region_prevalences
    )

  # Create accrual object
  accrual_instance <- accrual(
    treatment_arm_ids = trial_structure_instance@treatment_arm_ids,
    shared_control = shared_control,
    target_arm_size = target_arm_size,
    target_control = target_control,
    target_interim = target_interim,
    accrual_period = get_weeks(accrual_period),
    interim_period = get_weeks(interim_period),
    control_ratio = ctrl_ratio,
    centres_df = centres_df
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

  if (keep_files) {
    write.csv(
      accrual_instance@accrual, 
      paste0(output_path, "accrual-", run_time, ".csv"),
      row.names = FALSE
    )
  }

  # Plot outcome
  if (!quietly || keep_files) {
    p <- plot(accrual_instance)

    if (!quietly) {
      print(p)
    }
    if (keep_files) {
      ggplot2::ggsave(
        paste0(figs_path, "accrual-", run_time, ".png"),
        plot = p,
        width = 12,
        height = 8,
        dpi = 400
      )
    }
  }

  if (!quietly) {
    # Print accrual object
    print(accrual_instance)

    # Print summary of accrual object
    summary(accrual_instance)

    # Print trial structure object
    print(trial_structure_instance)

    cat("\n\nTreatment arm ids\n")
    print(trial_structure_instance@treatment_arm_ids_start)

    # Print summary of trial structure object
    cat("\n\nSite prevalences\n")
    summary(trial_structure_instance)$site_prev
  }

  # Plot trial structure object
  if (!quietly || keep_files) {
    p <- plot(trial_structure_instance)

    if (!quietly) {
      print(p)
    }

    if (keep_files) {
      ggplot2::ggsave(
        paste0(figs_path, "structure-", run_time, ".png"),
        plot = p,
        width = 12,
        height = 8,
        dpi = 400
      )
    }
  }

  # Return accrual object
  return(accrual_instance)
}
