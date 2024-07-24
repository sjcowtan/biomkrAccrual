
#' Driver for the procedure
#' @param target_arm_size Number of patients required per treatment arm
#' @param site_rates Vector of expected site accrual per month, one per site.
#' If all site rates are equal, one value can be given, if no_centres is
#' also specified.
#' @param no_centres Number of centres (defaults to the length of site_rates)
#' @param centre_starts Vector of site opening months, one per site.
#' if not specified, defaults to all starting at time zero.
#' @param fixed_site_rates Logical; are site rates to be considered as
#' fixed (TRUE) or do they vary randomly according to the gamma distribuiton
#' (FALSE); defaults to TRUE.
#' @param fixed_centre_starts Logical; do sites open at the beginning of 
#' the month specified (TRUE) or does the time vary randomly; defaults to TRUE.
#' @param site_caps Vector of site caps; if all sites are capped to the same
#' value can be a scalar. Defaults to NULL, for no site capping.
#' 
#' @examples 
#' sensitivity(10, 100)
#' 
#' @export
#' 
#' @importFrom rlang abort
#' 
sens_analysis <- function(
  target_arm_size = 308, 
  site_rates,
  no_centres = NULL,
  centre_starts = NULL,
  fixed_site_rates = TRUE,
  fixed_centre_starts = TRUE,
  site_caps = NULL,
  figs_path = "output_data/figures/"
) {

  # Input cleaning
  if (
    !is.numeric(target_arm_size) || 
      length(target_arm_size) != 1 ||
      (is.numeric(target_arm_size) && !identical(target_arm_size %% 1, 0))
  ) {
    rlang::abort("Target arm size should be one integer.")
  }

  if (!is.numeric(site_rates)) {
    rlang::abort("Site rates must be numeric.")
  }

  if (!is.null(no_centres) && !(is.numeric(no_centres))) {
    rlang::abort("If no_centres is specified, it must be an integer")
  } else if (!is.null(no_centres)) {
    if (!identical(no_centres %% 1, 0)) {
      rlang::abort("If no_centres is specified, it must be an integer")
    } else {
      no_centres <- as.integer(no_centres)
    }
  }

  if (length(site_rates) > 0) {
    if (is.null(no_centres)) {
      no_centres <- length(site_rates)
    } else if (length(site_rates) != no_centres && length(site_rates) != 1) {
      rlang::abort(paste(
        "Inconsistent number of centres; no_centres should", 
        "match the length of site_rates or be unspecified when",
        "more than one site rate is specified."
      ))
    }
  }

  if (
    !is.null(centre_starts) && (
      !is.numeric(centre_starts) || 
        (is.numeric(centre_starts) && !identical(centre_starts %% 1, 0))
    )
  ) {
    rlang::abort("Centre start times should be integers.")
  }

  # Choose appropriate analysis

  if (is.null(site_caps)) {
    # Do not need simulations
    if (fixed_centre_starts && fixed_site_rates && is.null(centre_starts)) {
      # Simplest case
      do_poisson_sensitivity(target_arm_size, site_rates, no_centres, figs_path)
    } else if (!fixed_centre_starts && is.null(centre_starts)) {
      do_poisson_gamma_sensitivity(
        target_arm_size, site_rates, no_centres, figs_path
      )
    }
  }
}


#' Sensitivity analysis for simplest case; fixed site rates, simultaneous centre
#' starts and no site capping.
#' 
#' @param target_arm_size Number of patients to be recruited.
#' @param site_rates Vector of site rates, one per site; or scalar if all 
#' identical.
#' @param no_centres Number of sites.
#' 
do_poisson_sensitivity <- function(target_arm_size, site_rates, no_centres, figs_path) {
  # Fixed, simultaneous centre starts and fixed site rates
  if (
    is.null(no_centres) || 
      identical(as.integer(no_centres), length(site_rates))
  ) {
    site_rates <- sum(site_rates)
  } else if (
    !is.null(no_centres) && identical(length(site_rates), as.integer(1))
  ) {
    site_rates <- site_rates * no_centres
  }
  
  do_sensitivity_plot_simultaneous( 
    target_arm_size, 
    site_rates,
    "sensitivity_poisson_",
    figs_path
  )

}

#' Sensitivity analysis for simple Poisson-gamma model; 
#' gamma-varying site rates, simultaneous centre starts 
#' and no site capping.
#' 
#' @param target_arm_size Number of patients to be recruited.
#' @param site_rates Vector of site rates, one per site; or scalar if all 
#' identical.
#' @param no_centres Number of sites.
#' 
#' @importFrom  stats rgamma
do_poisson_gamma_sensitivity <- function(
  target_arm_size, site_rates, no_centres, figs_path
) {

  site_rates <- sum(stats::rgamma(no_centres, rate = 1, shape = site_rates))

  do_sensitivity_plot_simultaneous(
    target_arm_size, 
    site_rates, 
    "sensitivity_poisson_gamma", 
    figs_path
  )

}


