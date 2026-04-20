#' @include prevalence.R



#' @title Command line 
#' 
#' @param shared_control TRUE if all treatment arms share the
#' same control arm; FALSE if each treatment arm has its own 
#' control. Defaults to TRUE.
#' @param target_times Vector of timings of interim and final
#' recruitment assessments, in months.
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
#' @param target_file Name of CSV file with target recruitment for each
#' arm at each interim analysis time and at final recruitment. This 
#' should have a column "arm" with the names of the treatment arms; 
#' this must match
#' the arm name as specified in in `arms_file`. Control targets are
#' not included as they can be deduced from the control_ratio and 
#' shared_control arguments. It may then have one or more columns for 
#' interim targets, and must have a "final" column for the final 
#' recruitment target.
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
#' @param seed Seed for random number generation.
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
  shared_control = TRUE,
  target_times = c(6, 12),
  precision = 10,
  var_lambda = 0.25,
  # active : control ratio (all active the same)
  control_ratio = c(1, 1),
  centres_file = "centres.csv",
  prop_file = "proportions.csv",
  target_file = "targets.csv",
  arms_file = "arms.json",
  data_path = "extdata/",
  output_path = "../biomkrAccrual_output_data/",
  figs_path = paste0(output_path, "/figures/"),
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE,
  fixed_region_prevalences = FALSE,
  quietly = FALSE,
  keep_files = TRUE,
  seed = 123456
) {
  
  checkmate::assert_logical(
    fixed_region_prevalences,
    any.missing = FALSE,
    null.ok = FALSE
  )

  checkmate::assert_logical(
    fixed_site_rates,
    any.missing = FALSE,
    null.ok = FALSE
  )

  checkmate::assert_logical(
    shared_control,
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

  checkmate::assert_numeric(
    control_ratio,
    any.missing = FALSE,
    lower = 10^-6,
    finite = TRUE,
    len = 2,
    null.ok = FALSE
  )

  if (!fixed_region_prevalences && is.null(precision)) {
    rlang::abort(paste("Either fixed_region_prevalences", 
      "must be TRUE, or a value must be given for the",
      "precision of the Dirichlet distribution."
    ))
  }


  if (fixed_site_rates && 
    checkmate::test_numeric(
      var_lambda, 
      any.missing = FALSE, 
      null.ok = FALSE
    )
  )  {
    rlang::warn(paste("Value given for var_lambda when", 
      "fixed_site_rates is TRUE: fixed_site_rates",
      "will take precendence."
    ))
    var_lambda <- NULL
  }

  checkmate::assert_numeric(
    target_times,
    any.missing = FALSE,
    lower = 0,
    finite = TRUE,
    min.len = 1,
    null.ok = FALSE
  )

  checkmate::assert_numeric(
    var_lambda,
    any.missing = FALSE,
    lower = 10^-6,
    finite = TRUE,
    len = 1,
    null.ok = TRUE
  )

  if (!fixed_site_rates && is.null(var_lambda)) {
    rlang::abort(paste("Either fixed_site_rates", 
      "must be TRUE, or a value must be given for the",
      "variance of the site rates."
    ))
  }

  checkmate::assert_integerish(
    seed,
    min.len = 1,
    max.len = 7,
    null.ok = TRUE
  )


  # Verify inputs
  ## append "/" if no slash on end
  data_path <- gsub("(\\w+)$", "\\1/", data_path)
  output_path <- gsub("(\\w+)$", "\\1/", output_path)
  figs_path <- gsub("(\\w+)$", "\\1/", figs_path)


  # Make into full path so only one set of syntax needed
  if (!grepl("^/", output_path)) {
    output_path <- paste0(getwd(), "/", output_path)
  }
  
  ## Check data directory exists
  if (grepl("^extdata/?$", data_path)) {
    checkmate::assert_directory_exists(system.file(
      data_path, package = "biomkrAccrual"
    ), access = "rx")
  } else {
    checkmate::assert_directory_exists(data_path, access = "rx")
  }

  # Set up output directory if does not already exist
  makeifnot_dir(output_path)
  
  # Set up output figures directory if does not already exist
  makeifnot_dir(figs_path)

  # Switch between system.file() if using data shipped with the
  # the package, or straight file access if now.

  if (grepl("^extdata/?$", data_path)) {
    prop_file <- system.file(
      data_path, prop_file, package = "biomkrAccrual"
    )
    arms_file <- system.file(
      data_path, arms_file, package = "biomkrAccrual"
    )
    centres_file <- system.file(
      data_path, centres_file, package = "biomkrAccrual"
    )
    target_file <- system.file(
      data_path, target_file, package = "biomkrAccrual"
    )
  } else {
    prop_file <- paste(data_path, prop_file, sep = "/")
    arms_file <- paste(data_path, arms_file, sep = "/")
    centres_file <- paste(data_path, centres_file, sep = "/")
    target_file <- paste(data_path, target_file, sep = "/")
  }

  # Check input files exist and are readable
  checkmate::assert_file_exists(prop_file, access = "r")
  checkmate::assert_file_exists(arms_file, access = "r")
  checkmate::assert_file_exists(centres_file, access = "r")
  checkmate::assert_file_exists(target_file, access = "r")

  
  # Read parameters
  prop_params_df <- utils::read.csv(prop_file) 
    
  arms_ls <- 
    jsonlite::read_json(arms_file, simplifyVector = TRUE)
  
  centres_df <- utils::read.csv(centres_file)

  target_df <- utils::read.csv(target_file)

  # Add fail if read fails

  # Fail if centres_file is in wrong format
  if (!(
    isTRUE(all.equal(
      names(centres_df), 
      c("site", "start_month", "mean_rate", "region", "site_cap")
      # site_cap is optional, no cap if not present
    )) || isTRUE(all.equal(
      names(centres_df), c("site", "start_month", "mean_rate", "region")
    ))
  )) {
    rlang::abort(paste(
      "Format error: centres.csv should have columns site,",
      "start_month, mean_rate, region, and optionally site_cap"
    ))
  }

  # Fail if arms file is in wrong format
  if (any(duplicated(names(arms_ls)))) {
    rlang::abort(paste(
      "Format error: arm names duplicated in arms file"
    ))
  }

  # Fail if target_file is in wrong format
  if (!(
    isTRUE(all(
      c("arm", "final") %in% names(target_df)
    ))
  )) {
    rlang::abort(
      "Format error: target.csv must have columns target and final"
    )
  }
  if (ncol(target_df) - 1 != length(target_times)) {
    rlang::abort(paste(
      "Format error: target_times must have one time value for each",
      "target column in the targets file" 
    ))
  }
  if (isFALSE(checkmate::testSetEqual(names(arms_ls), target_df$arm))) {
    # Tests for names matching in any order
    rlang::abort(paste(
      "Format error: arm names in the targets file are not the same",
      "as those in the arms file"
    ))
  }
  if (any(duplicated(target_df$arm))) {
    rlang::abort(paste(
      "Format error: arm names duplicated in targets file"
    ))
  }

  # Set run_time to timestamp output files
  run_time <- format(Sys.time(), "%F-%H-%M-%S")

  # Set seed
  if (is.null(seed)) {
    if (length(.Random.seed) != 7 || .Random.seed[1] != 10407) {
      rlang::abort("No l'Ecuyer seed set and seed is NULL.")
    }
  } else if (length(seed) == 1) {
    set.seed(seed)
  } else if (length(seed) == 7 && seed[1] == 10407) {
    # Set seed for next L'Ecuyer stream
    ##### think about this!  Have not set kind
    assign(
      x = ".Random.seed", 
      # This still isn't making the output different,
      # but there is variation where there shouldn't be
      # Seeds not propagating properly
      value = seed,
      envir = as.environment(-1)
    )
  } else {
    rlang::abort("Invalid seed, must be L'Ecuyer or a single integer.")
  }

  # Get start weeks & order centres_df by start week and site number
  centres_df <- do_clean_centres(centres_df)
  centres_df$start_week <- get_weeks(centres_df$start_month - 1) + 1

  # Sort target times and convert to weeks
  target_times <- get_weeks(sort(target_times))

  # Sort target_df arm columns by value with the arm name first
  arm_col <- which(colnames(target_df) == "arm")

  target_df <- 
    target_df[, c(1, 1 + order(unlist(target_df[1, -1])))]

  # Sort target_df by rows to match order of arms_ls
  if (any(target_df$arm != names(arms_ls))) {
    target_df <- target_df[match(names(arms_ls), target_df$arm), ]
  }

  # Total target recruitment
  target_recruit <- round(
    sum(control_ratio) * sum(target_df$final) / 
      control_ratio[1], 0
  )

  # Make control ratio sum to 1
  control_ratio <- control_ratio / sum(control_ratio)

  # Complete site cap if incomplete, using total recruitment target
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
      control_ratio = control_ratio,
      fixed_region_prevalences = fixed_region_prevalences
    )

  # Create accrual object
  accrual_instance <- accrual(
    treatment_arm_ids = trial_structure_instance@treatment_arm_ids,
    shared_control = shared_control,
    target_df = target_df,
    target_times = target_times,
    control_ratio = control_ratio,
    fixed_site_rates = fixed_site_rates,
    var_lambda = var_lambda,
    centres_df = centres_df
  )

  while (
    # Any arms are recruiting
    any(trial_structure_instance@treatment_arm_struct) &&
      # Any sites are recruiting
      length(accrual_instance@active_sites) > 0 
  ) {

    # Add a week's accrual
    obj_list <- accrue_week(
      accrual_instance, 
      trial_structure_instance
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
