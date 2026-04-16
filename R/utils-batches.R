#' Set seeds for multiple simulations, using a new L'Ecuyer 
#' stream for each one.
#' 
#' @param n No. runs
#' @param seed Starting seed
#' @return List of N sets of 7 seeds
#' @importFrom parallel nextRNGStream
#' 

getseeds <- function(
  n = 5, 
  seed = NULL
) {

  # Get current RNG and seed for current environment
  oldrng <- RNGkind()[1L]
  oldseed <- get(x = ".Random.seed", envir = as.environment(-1))

  # On exit: Reset RNG and seed for current environment to previous value
  on.exit(
    set.seed(
      kind = oldrng, 
      seed = oldseed
    )
  )
  if (is.null(seed)) {
    seed <- oldseed
  } 

  # Make sure the seed is at least plausible
  stopifnot(
    is.numeric(seed), 
    all(is.finite(seed))
  )

  # Set a L'Ecuyer seed
  set.seed(
    kind = "L'Ecuyer-CMRG", 
    seed = seed
  )
  le_seed <- get(".Random.seed", envir = as.environment(-1))

  # Initialise empty seed list, one per simulation
  seeds <- vector(mode = "list", length = n)

  for (i in 1:n) {
    le_seed <- parallel::nextRNGStream(le_seed)
    seeds[[i]] <- le_seed
  }

  return(seeds)
}


#' Run a number of batches of recruitment prediction and
#' collect summary statistics on arm closures, and final
#' recruitment totals for all experimental and control arms.
#' Seeds for each run are generated using the L'Ecuyer 
#' pseudo-random number generator.
#' 
#' @param n Number of instances to run (defaults to 100)
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
#' should have a column "target" with the arm name; this must match
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
#' @param seed Seed to be used to generate the seeds for each simulation.
#' @return Dataframe of site closing times
#' @return Dataframe of experimental arm totals
#' @export
#' 
#' @importFrom jsonlite read_json
#' @importFrom utils write.csv
#' 
biomkrAccrualSim <- function(
  n = 100,
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
  figs_path = paste0(output_path, "figures/"),
  fixed_centre_starts = TRUE,
  fixed_site_rates = FALSE,
  fixed_region_prevalences = FALSE,
  quietly = FALSE,
  keep_files = TRUE,
  seed = 123456
) {
  # Timestamp for batch files (but not individual run files)
  run_time <- format(Sys.time(), "%F-%H-%M-%S")

  ## Check data directory exists
  if (grepl("^extdata/?$", data_path)) {
    checkmate::assert_directory_exists(system.file(
      data_path, package = "biomkrAccrual"
    ), access = "rx")
  } else {
    checkmate::assert_directory_exists(data_path, access = "rx")
  }

  if (grepl("^extdata/?$", data_path)) {
    jsonfile <- system.file(
      data_path, arms_file, package = "biomkrAccrual"
    )
  } else {
    jsonfile <- paste(data_path, arms_file, sep = "/")
  }

  # Information for setting up dataframe
  arms_ls <- 
    jsonlite::read_json(jsonfile, simplifyVector = TRUE)

  no_arms <- ifelse(
    shared_control,
    length(arms_ls) + 1,
    2 * length(arms_ls)
  )

  # Define matrix of zeroes for efficiency
  arm_closures_mx <- structure(
    matrix(0, nrow = n, ncol = ifelse(
      shared_control,
      length(arms_ls) + 1,
      2 * length(arms_ls)
    )),
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
  # List with one element for each assessment time
  arm_interim_ls <- lapply(
    seq_len(length(target_times)),
    function(x) {
      structure(
        matrix(
          0, 
          nrow = n,
          ncol = ifelse(
            shared_control,
            length(arms_ls) + 1,
            2 * length(arms_ls)
          )
        ),
        class = c("armtotals", "matrix", "array")
      )
    }
  )
  arm_accrual_ls <- list(length = n)

  # Set column names
  if (shared_control) {
    colnames(arm_totals_mx) <- c(names(arms_ls), "Control")
  } else {
    colnames(arm_totals_mx) <- c(names(arms_ls), paste0("C-", names(arms_ls)))
  }
  colnames(arm_closures_mx) <- colnames(arm_totals_mx)
  arm_interim_ls <- 
    lapply(arm_interim_ls, "colnames<-", colnames(arm_totals_mx))

  # Get current RNG and seed for current environment
  oldrng <- RNGkind()
  oldseed <- get(x = ".Random.seed", envir = as.environment(-1))

  # On exit: Reset RNG and seed for current environment to previous value
  on.exit(
    RNGkind(
      oldrng[1L], 
      normal.kind = oldrng[2L], 
      sample.kind = oldrng[3L]
    )
  )
  on.exit(
    set.seed(
      kind = oldrng[1L], 
      seed = oldseed
    ),
    add = TRUE
  )

  # Generate a list of seeds, one per simulation
  seeds <- getseeds(n, seed)

  # Set the random number generator to the one compatible
  # with RNGstreams
  RNGkind(
    kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", 
    sample.kind = "Rejection"
  )

  # Run batches
  for (irun in seq(n)) {

    # Set seed for next L'Ecuyer stream
    assign(
      x = ".Random.seed", 
      value = seeds[[irun]],
      envir = as.environment(-1),
      inherits = TRUE
    )

    # Run one simulation
    accrual_instance <- biomkrAccrual(
      shared_control = shared_control,
      target_times = target_times,
      precision = precision,
      var_lambda = var_lambda,
      control_ratio = control_ratio,
      centres_file = centres_file,
      prop_file = prop_file,
      target_file = target_file,
      arms_file = arms_file,
      data_path = data_path,
      fixed_centre_starts = fixed_centre_starts,
      fixed_site_rates = fixed_site_rates,
      fixed_region_prevalences = fixed_region_prevalences,
      quietly = TRUE,
      keep_files = FALSE,
      seed = NULL
    )

    # Tally arm closure times
    arm_closures_mx[irun, ] <- accrual_instance@phase_changes
    arm_totals_mx[irun, ] <- treat_sums(
      accrual_instance@accrual[
        seq_len(min(
          nrow(accrual_instance@accrual), 
          accrual_instance@target_times[length(accrual_instance@target_times)]
        )), ,
      ]

    )

    # Tally arm totals at each interim analysis & planned or actual end
    for (i in seq_len(length(accrual_instance@target_times))) {
      arm_interim_ls[[i]][irun, ] <- treat_sums(
        accrual_instance@accrual[
          seq(min(
            accrual_instance@target_times[i],
            dim(accrual_instance@accrual)[1]
          )), , 
        ]
      )
    }

    # Don't know how many weeks to predeclare => not array. List of matrices
    arm_accrual_ls[[irun]] <- apply(
      accrual_instance@accrual,
      1:2,
      sum
    )
  }

  # Reformat accrual list to be all the same length
  most_weeks <- max(unlist(lapply(arm_accrual_ls, nrow)))
  accrual_ls <- lapply(1:n, function(x) {
    matrix(
      0, 
      nrow = most_weeks,
      ncol = no_arms
    )
  })
  
  for (i in seq_len(n)) {
    accrual_ls[[i]][seq_len(nrow(arm_accrual_ls[[i]])), ] <-
      arm_accrual_ls[[i]]
  }

  # Now convert to array
  arm_accrual_ar <- simplify2array(accrual_ls)
  dimnames(arm_accrual_ar)[[2]] <- dimnames(accrual_instance@accrual)[[2]]
  names(dimnames(arm_accrual_ar)) <- c("Week", "Arm", "Simulation")


  # And now to a list of matrices by arm
  accrual_byarm_ls <- asplit(arm_accrual_ar, 2)
  accrual_byarm_ls <- lapply(
    accrual_byarm_ls,
    function(m) structure(m, class = c("armaccrual", "matrix", "array"))
  )

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

  for (i in seq_len(length(arm_interim_ls))) {
    write.csv(
      as.data.frame(arm_interim_ls[[i]]),
      paste0(
        output_path, "arm_interim_totals_", target_times[i], 
        "mo_", datetime, ".csv"
      ),
      row.names = FALSE
    )
  }

  # Write simulation data as a JSON of a list, one element per
  # arm, consisting of matrices of simulation number * week
  jsonlite::write_json(
    accrual_byarm_ls,
    paste0(output_path, "arm_accrual_", datetime, ".json")
  )

  print("Arm closure weeks")
  print(summary(as.data.frame(arm_closures_mx)))
  print("Arm totals")
  print(summary(as.data.frame(arm_totals_mx)))
  for (i in seq_len(length(target_times))) {
    print(paste0(target_times[i], "mo. accrual:"))
    print(summary(as.data.frame(arm_interim_ls[[i]])))
  }

  # Going to need to know recruitment targets for plots
  ### This has been read before so shouldn't fail
  if (grepl("^extdata/?$", data_path)) {
    target_file <- system.file(
      data_path, target_file, package = "biomkrAccrual"
    )
  } else {
    target_file <- paste(data_path, target_file, sep = "/")
  }

  target_df <- utils::read.csv(target_file)
  target_expanded_df <- expand_targets(target_df)

  # Accrual for all arms at each analysis point
  for (i in seq_len(length(arm_interim_ls))) {
    p <- plot(
      arm_interim_ls[[i]], 
      target = unique(target_expanded_df[, i + 1]), 
      target_names = target_group(target_expanded_df, target_col = i + 1),
      target_week = target_times[i],
      plot_id = paste("Accrual at", target_times[i], "months")
    )

    ggplot2::ggsave(
      paste0(
        figs_path, "arm-totals-", target_times[i], 
        "mo-", run_time, ".png"
      ),
      plot = p,
      width = 12,
      height = 8,
      dpi = 400
    )

    print(p)
  }

  p <- plot(
    arm_totals_mx, 
    target = unique(target_expanded_df[, ncol(target_expanded_df)]), 
    target_names = target_group(
      target_expanded_df,
      target_col = ncol(target_expanded_df)
    ),
    plot_id = "Total accrual"
  )

  ggplot2::ggsave(
    paste0(figs_path, "arm-totals-", run_time, ".png"),
    plot = p,
    width = 12,
    height = 8,
    dpi = 400
  )

  print(p)
  

  # Individual accrual plots
  arm_names <- dimnames(arm_totals_mx)[[2]]

  ## Mark arms as treatment or control
  treatment_arms <- startsWith(arm_names, "T")

  ## Same colours as in interim plot
  col_order <- c(seq_len(length(treatment_arms))[-1], 1)
  palette <- grDevices::palette.colors(
    palette = "R4",
    length(treatment_arms) + 1
    # One pink is more than enough
  )[-6]
  arm_colours <- palette[col_order]

  
  # Total accrual plots

  # Loop across all arms
  for (i in seq_len(length(accrual_byarm_ls))) {
    ## Loop across interim and total
    for (j in seq_len(ncol(target_expanded_df[, -1]))) {
    
      p <- accrual_arm_plot(
        arm_interim_ls[[j]],
        arm_colours, 
        treatment_arms,
        target_expanded_df[[i, j + 1]],
        plot_id = paste(target_times[j], "month"),
        i
      )

      ggplot2::ggsave(
        paste0(
          figs_path, 
          "arm-totals-",
          target_times[j],
          "mo-",
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

  ## Plot simulations over time for each arm

  for (i in seq_len(length(accrual_byarm_ls))) {
    target_index <- 2 - as.numeric(treatment_arms[[i]])
    name <- names(accrual_byarm_ls)[[i]]
    p <- plot(
      accrual_byarm_ls[[i]],
      arm_colour = arm_colours[i],
      target = unlist(target_expanded_df[i, -1]),
      target_names = paste(target_times, "month"),
      plot_id = name
    )
    ggplot2::ggsave(
      paste0(
        figs_path, 
        "arm-accrual-",
        tolower(names(accrual_byarm_ls)[j]),
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

  # Plot time to accrual for each target and each arm
  for (i in seq_len(length(accrual_byarm_ls))) {
    # Data extraction for interim and accrual
    accrual_times <- threshold_week(
      accrual_byarm_ls[[i]], 
      unlist(target_expanded_df[i, -1])
    )

    for (j in seq_len(length(target_times))) {
      p <- plot(
        accrual_times[[j]],
        arm_colour = arm_colours[i],
        target = unlist(target_expanded_df[i, j + 1]),
        target_names = paste(target_times[j], "month"),
        plot_id = names(accrual_byarm_ls)[[i]]
      )
      ggplot2::ggsave(
        paste0(
          figs_path, 
          "arm-week-", 
          names(accrual_byarm_ls)[[i]],
          "-",
          target_times[j], 
          "mo-",
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

