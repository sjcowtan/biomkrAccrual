#' A class to contain the accrued patient allocations,
#' by site, treatment/control arm and week
#' 
#' @slot accrual 3-D array with axes site, experimental arm
#' and week
#' @slot target_df Dataframe containing the targets for the number of 
#' patients recruited by each arm at each timepoint in `target_times`.
#' @slot target_times Vector of times of analyses or other recruitment
#' progress assessments (weeks).
#' @slot phase_changes Vector of week numbers when arms closed
#' @slot site_closures Vector of weeks sites closed; NA indicates open
#' @slot week Current recruitment week
#' @slot active_arms Vector of indices of open arms
#' @slot active_sites Vector of indices of open sites (automatically
#' generated)
#' @slot closed_sites Vector of indices of closed sites
#' @slot shared_control TRUE if a shared control arm is being used, else FALSE
#' @slot fixed_site_rates TRUE if expected site rate to be used; FALSE 
#' draws the site rate from a gamma distribution
#' @slot site_in_region Vector of indices for which set of expected 
#' prevalences each site should use 
#' @slot site_cap Vector of maximum number of patients for each site
#' @slot site_mean_rate Vector of expected recruitment rates for each site
#' @slot site_rate Vector of recruitment rates for each site, drawn from a 
#' gamma distribution
#' @slot start_week Vector of weeks that each site opens recruitment
#' @slot site_index Vector of index numbers for each site
#' @slot treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @slot var_lambda Variance of site recruitment rates.
#' 
#' @param treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @param shared_control TRUE if all experimental arms share one control arm;
#' FALSE if each has their own
#' @param fixed_site_rates TRUE if expected site rate to be used; FALSE 
#' draws the site rate from a gamma distribution
#' @param target_df Dataframe containing the targets for the number of 
#' patients recruited by each arm at each timepoint in `target_times`.
#' @param target_times Vector of times of analyses or other recruitment
#' progress assessments (weeks).
#' @param control_ratio Ratio of patient allocation to treatment arm
#' versus control for all active arms; defaults to c(1, 1).
#' @param var_lambda Variance of site recruitment rates.
#' @param centres_df Dataframe with columns "site", "start_month", "mean_rate", 
#' "region" and "site_cap"
#' 
#' @usage accrual(treatment_arm_ids, shared_control, fixed_site_rates, 
#' target_df, target_times, control_ratio, var_lambda, centres_df)
#' @name accrual
#'
#' @export
#' 
accrual <- S7::new_class("accrual",
  package = "biomkrAccrual",
  properties = list(
    accrual = S7::class_integer,
    target_times = S7::class_integer,
    target_df = S7::class_data.frame,
    control_ratio = S7::class_double,
    phase_changes = S7::class_integer,
    site_closures = S7::class_integer,
    week = S7::class_integer,
    active_arms = S7::class_integer,
    shared_control = S7::class_logical,
    fixed_site_rates = S7::class_logical,
    site_in_region = S7::class_integer,
    site_cap = S7::class_integer,
    site_mean_rate = S7::class_double,
    site_rate = S7::class_double,
    site_start_week = S7::class_integer,
    site_index = S7::class_integer,
    treatment_arm_ids = S7::class_list,
    closed_sites = S7::class_integer,
    active_sites = S7::new_property(
      getter = function(self) {
        setdiff(
          which(self@site_start_week <= self@week),
          self@closed_sites
        )
      }
    )
  ),
  constructor = function(
    treatment_arm_ids = S7::class_missing,
    shared_control = S7::class_missing,
    fixed_site_rates = S7::class_missing,
    target_df = S7::class_missing,
    target_times = S7::class_missing,
    control_ratio = S7::class_missing,
    var_lambda = S7::class_missing,
    centres_df = S7::class_missing
  ) {
    # Create the object and populate it
    S7::new_object(
      # Parent class
      S7::S7_object(),
      # weeks * treatments * no_centres
      accrual = array(
        # Defaults to 0 for no patients recruited
        as.integer(0), 
        dim = c(
          # Max. weeks
          target_times[length(target_times)],
          # No. experimental arms including control
          length(treatment_arm_ids) + 
            ifelse(shared_control, 1, length(treatment_arm_ids)),
          # No. sites 
          length(unique(centres_df$site))
        ),
        dimnames = list(
          Weeks = NULL,
          Arms = c(
            names(treatment_arm_ids), 
            if (shared_control) {
              "Control"
            } else {
              paste("Control", names(treatment_arm_ids))
            }
          ),
          Centres = c(paste("Centre", unique(centres_df$site)))
        )
      ),
      shared_control = shared_control,
      control_ratio = control_ratio,
      target_df = expand_targets(target_df, shared_control, control_ratio),
      target_times = as.integer(target_times),
      treatment_arm_ids = treatment_arm_ids,
      phase_changes = rep(
        NA_integer_, 
        length(treatment_arm_ids) + 
          ifelse(shared_control, 1, length(treatment_arm_ids))
      ),
      site_closures = rep(NA_integer_, length(unique(centres_df$site))),
      week = as.integer(1),
      active_arms = seq(length(treatment_arm_ids) +
        ifelse(
          shared_control,
          1,
          length(treatment_arm_ids)
        )
      ),
      fixed_site_rates = fixed_site_rates,
      site_in_region = as.integer(centres_df$region),
      site_cap = as.integer(centres_df$site_cap),
      site_mean_rate = centres_df$mean_rate,
      site_start_week = as.integer(centres_df$start_week),
      site_index = as.integer(centres_df$site),
      site_rate = set_site_rates(
        mean_rate = centres_df$mean_rate, 
        fixed_site_rates = fixed_site_rates, 
        var_lambda = var_lambda
      ),
      closed_sites = NA_integer_
    )
  }
)


#' Sum accrual array by site, to enable checking against site caps.
#' @param obj Object of class "accrual".
#' @param ... For R CMD check compatibility.
#' 
#' @return vector of total accrual by recruitment site.
#' 
#' @export 
#' 
site_sums <- S7::new_generic("site_sums", "obj")
S7::method(site_sums, accrual) <- function(obj) {
  # Sum across dimensions 1:dims
  colSums(obj@accrual, dims = 2)
}


# Sum accrual array by experimental arm (including control).
#' 
#' @param x Accrual array with dimensions Weeks, Arms and Centres.
#' @param control_total Logical; if TRUE return single total for all 
#' control arms (not used if `shared_control` is TRUE); defaults to FALSE.
#' @param ... When called on an array the arguments no_treat_arms and 
#' shared_control may be necessary. Additional arguments are 
#' unnecessary if calling treat_sums.biomkrAccrual::accrual().
#' @export
treat_sums <- function(x, control_total, ...) {
  UseMethod("treat_sums", x)
}


#' Sum accrual array by experimental arm (including control).
#' 
#' @param x Accrual array with dimensions Weeks, Arms and Centres.
#' @param control_total Logical; if TRUE return single total for all 
#' control arms (not used if `shared_control` is TRUE); defaults to FALSE.
#' @param ... Additional arguments (none needed).
#' @param no_treat_arms Number of treatment arms (as opposed to control 
#' arms); only needed if shared_control is FALSE. MUST use 
#' `no_treat_sums =` to set.
#' @param shared_control TRUE if all treatment arms share the
#' same control arm; FALSE if each treatment arm has its own 
#' control. Defaults to TRUE. MUST use `shared_control =` to set.
#' 
#' @return vector of total accrual by experimental arm.
#' 
#' @importFrom checkmate assert_logical assert_integerish
#' @export
#' 
treat_sums.array <- function(
  x, 
  control_total = FALSE,
  ...,
  no_treat_arms,
  shared_control = TRUE
) {

  # Check that arguments are consistent
  if (!shared_control) {
    checkmate::assert_logical(control_total)
    checkmate::assert_integerish(
      no_treat_arms,
      lower = 1,
      any.missing = FALSE,
      len = 1,
      null.ok = FALSE
    )
  }

  # Permute the array so that the first dimension is the
  # dimension you want to get sums for (experimental arms)
  arm_sums <-
    as.integer(colSums(rowSums(x, dims = 2, na.rm = TRUE)))

  # If want total for control arms rather than separate values
  if (!shared_control && control_total) {
    arm_sums <- c(
      # Treatment arms
      arm_sums[seq_along(no_treat_arms)],
      # Control
      sum(arm_sums[seq(no_treat_arms + 1, length.out = no_treat_arms)])
    )
  }

  return(arm_sums)
}


#' Sum accrual array element of accrual object, by experimental 
#' arm (including control).
#' 
#' @name treat_sums
#' 
#' @param x Object of class `accrual`. 
#' @param control_total Logical; if TRUE return single total for all 
#' control arms
#' @param ... Additional arguments (none needed).
#' 
#' @return vector of total accrual by experimental arm
#' 
#' @export
#' 
`treat_sums.biomkrAccrual::accrual` <- function(
  x,
  control_total = FALSE,
  ...
) {

  # Call treat_sums.array() on accrual array element
  treat_sums(
    x@accrual,
    control_total, 
    no_treat_arms = length(x@phase_changes), 
    shared_control = x@shared_control,
  )
}


#' Select arms to cap for single site, or sites to cap for single arm
#' @param population Vector of allocations of patients this week
#' @param captotal Value of cap for site or arm (as appropriate)
#' @return Vector of arms to cap (can include multiples of the same arm)
#' 
do_choose_cap <- function(population, captotal) {
  if (length(population) == captotal) { 
    # Use whole population if correct length
    capped <- population
  } else {  
    # Sample    
    capped <- sample(
      population, 
      size = min(length(population), captotal),
      replace = FALSE
    )
  }

  return(capped)
}


#' Get no. weeks from months
#' Currently 4 week month, fix later using lubridate
#' @param months Duration in months
#' @return Duration in weeks
#' 
#' @importFrom checkmate assert_numeric
#' 
get_weeks <- function(months) {

  checkmate::assert_numeric(
    months,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1,
    null.ok = FALSE
  )
  as.integer(round(months * 4, 0))
}


#' Generate site rates for each site, drawing from a gamma 
#' distribution if the rates are not fixed. 
#' 
#' @param mean_rate Vector of expected recruitment rates for all sites,
#' in patients per month (from the centres file). 
#' @param fixed_site_rates Logical: whether site rates are fixed or 
#' distributed according to the gamma distribution.
#' @param var_lambda Variance of the gamma distribution; NULL if 
#' site rates are fixed.
#' 
#' @return Vector of site rates in patients per week.
#' 
#' @importFrom stats rgamma
#' 
set_site_rates <- function(mean_rate, fixed_site_rates, var_lambda) {

  # mean_rate is in recruitment per month, convert to weeks
  if (fixed_site_rates) {
    rates <- mean_rate / get_weeks(1)
  } else {
    rates <- stats::rgamma(
      n = length(mean_rate),
      shape = mean_rate^2 / var_lambda,
      # Per week not per month
      rate = mean_rate / var_lambda
    ) / get_weeks(1)
  }

  return(rates)
}


#' Implement site cap on a week's accrual. 
#' @param obj Accrual object
#' @param ... For compliance with R CMD check.
#' 
#' @return Modified accrual object with capped week's accrual and 
#' with any capped sites removed from active_sites
#' 
apply_site_cap <- S7::new_generic("apply_site_cap", "obj")
S7::method(apply_site_cap, accrual) <- function(obj, ...) {
  # If any sites exceed their cap, remove accrual from randomly
  # selected arms until sites are at cap 
  # Represents sites closing during the week
  site_captotal <- site_sums(obj) - obj@site_cap
  if (any(site_captotal > 0)) {
    # Loop over sites which are above the cap
    for (site in which(site_captotal > 0)) {
      # Vector of instances of populated arms in week's accrual, 
      # including control, e.g. c(1, 1, 2, 4)
      population <- unlist(sapply(
        which(obj@accrual[obj@week, , site] > 0), 
        function(j) rep(j, obj@accrual[obj@week, j, site])
      ))
      
      # Randomly select population instances to remove, leaving
      # enough to max out the cap
      if (length(population) > 0) {
        capped <- do_choose_cap(population, site_captotal[site])

        # Remove those instances
        for (i in capped) {
          obj@accrual[obj@week, i, site] <- 
            as.integer(obj@accrual[obj@week, i, site] - 1)
        }
      }
    }
  }

  # Don't remove sites from active_sites here, because they may
  # be reinstated during arm capping
  
  return(obj)
}


#' Implement arm cap on week's accrual to experimental arms
#' @param accrual_obj Object of class `accrual`
#' @param struct_obj Object of class `trial_structure`
#' @param ... For compliance with R CMD check.
#' @return Modified accrual object with capped week's accrual
#' 
#' @importFrom rlang abort
#' 
apply_arm_cap <- 
  S7::new_generic("apply_arm_cap", c("accrual_obj", "struct_obj"))
S7::method(apply_arm_cap, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj, ...) {
    
    # Get totals for experimental arms (dropping control)
    arm_sums <- 
      #treat_sums(accrual_obj)[seq_len(length(accrual_obj@phase_changes))]
      treat_sums(accrual_obj)
  
    # Compare with cap
    arm_captotal <- arm_sums - accrual_obj@target_df$final
    over_cap_all <- arm_captotal >= 0
    over_cap <- over_cap_all[accrual_obj@active_arms]

    active_tocap <- NULL
    
    if (sum(over_cap) > 0) {
      if (accrual_obj@shared_control) {
        # If any arms remain, there must be control and
        # at least one experimental arm
        if (xor(
          # Any active arms remaining open if these are closed
          any(
            !over_cap[-length(over_cap)],
            na.rm = TRUE
          ),
          # Shared control being closed
          over_cap[length(over_cap)]
        )) {
          # Safe to close all at cap
          active_tocap <- accrual_obj@active_arms[over_cap]
        } else if (
          # At least one active arm to close
          sum(over_cap) >= 2 
          &&
            # Some active arms would remain open
            any(
              !over_cap[-length(over_cap)],
              na.rm = TRUE
            )
          && 
            # Control at or over cap
            over_cap[length(over_cap)]
        ) {
          # Close the active arms but not control
          active_tocap <- 
            accrual_obj@active_arms[over_cap][
              -length(accrual_obj@active_arms[over_cap])
            ]
        } 
      } else {
        active_tocap <- 
          which(colSums(matrix(
            over_cap_all, 
            nrow = 2, 
            byrow = TRUE
          )) == 2)
        # Cap control arms as well as active
        active_tocap <- c(
          active_tocap, 
          active_tocap + (length(over_cap_all) / 2)
        )
      }
    }

    # Apply cap to affected arms
    if (!is.null(active_tocap)) {
      for (arm in active_tocap) {
        # Which sites recruited to that arm this week?
        population <- as.integer(unlist(sapply(
          which(accrual_obj@accrual[accrual_obj@week, arm, ] > 0), 
          function(j) rep(j, accrual_obj@accrual[accrual_obj@week, arm, j])
        )))

        # Possible accruals to remove
        if (length(population) > 0) {
          capped <- do_choose_cap(population, arm_captotal[arm])
          # Remove capped instances
          for (i in capped) {
            accrual_obj@accrual[accrual_obj@week, arm, i] <- 
              as.integer(accrual_obj@accrual[accrual_obj@week, arm, i] - 1)
          }
        }
      }
    } 

    # Record closing week for capped arms
    accrual_obj@phase_changes[active_tocap] <- 
      accrual_obj@week
    
    # Update active_arms
    accrual_obj@active_arms <- setdiff(accrual_obj@active_arms, active_tocap)

    # Also on the trial structure object
    struct_obj <- remove_treat_arms(
      struct_obj, 
      arms = active_tocap[active_tocap <= ncol(struct_obj@treatment_arm_struct)]
    )

    # Recheck site caps
    site_captotal <- site_sums(accrual_obj) - accrual_obj@site_cap
    
    # Record closing week for capped sites
    capped_sites <- site_captotal[accrual_obj@active_sites] >= 0
    accrual_obj@site_closures[accrual_obj@active_sites[capped_sites]] <-
      accrual_obj@week

    # Update active sites(accrual_obj)
    accrual_obj@closed_sites <- which(site_captotal >= 0)

    return(list(accrual_obj, struct_obj)) 
  }

#' Randomise a week's expected accrual amongst the sites, according to 
#' prevalence.
#' @param accrual_obj An object of class `accrual`.
#' @param struct_obj An object of class `trial_structure`.
#' be treated as exact; FALSE if they should be drawn from a gamma
#' distribution with a mean of the specified rate.
#' @param ... For R CMD check compatibility.
#' 
#' @return Matrix of week's accrual by site and recruitment arm.
#' 
week_accrue <- S7::new_generic("week_accrue", c("accrual_obj", "struct_obj"))
S7::method(week_accrue, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj) {

    # Update the site rates
    #accrual_obj <- set_site_rates(accrual_obj)

    # Initialising (sites * experimental arms)
    week_mx <- matrix(
      as.integer(0), 
      nrow = dim(accrual_obj@accrual)[2], 
      ncol = dim(accrual_obj@accrual)[3]
    )

    # Accrual per site is poisson distributed

    ## Get prevalences for active arms by site
    prev_mx <- colSums(
      struct_obj@experimental_arm_prevalence[, , accrual_obj@active_sites]
    )

    if (is.matrix(prev_mx)) {
      # multiply prevalences by lambdas for each site
      # => site * arm matrix
      lambda_prev_mx <- sweep(
        prev_mx, 
        MARGIN = 2, 
        accrual_obj@site_rate[accrual_obj@active_sites], 
        `*`
      )
    } else {
      lambda_prev_mx <- accrual_obj@site_rate[accrual_obj@active_sites] *
        prev_mx
    }

    ## Draw from Poisson distribution
    week_mx[, accrual_obj@active_sites] <-
      matrix(
        rpois(
          n = length(lambda_prev_mx), 
          lambda = lambda_prev_mx
        ), 
        ncol = length(accrual_obj@active_sites)
      )

    return(list(accrual_obj, week_mx))
  }


#' Generate accrual for a week and assign it to the 
#' accrual object. Increments the property "week",
#' which holds the next week number to accrue.
#' @param accrual_obj An object of class "accrual"
#' @param struct_obj An object of class "trial_structure"
#' @param ... For R CMD check compatibility.
#' 
#' @return An object of class "accrual"
#'
accrue_week <- S7::new_generic("accrue_week", c("accrual_obj", "struct_obj"))
S7::method(accrue_week, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj, ...) {

    # Should not get here if there aren't any but
    if (length(accrual_obj@active_sites) > 0) {

      week_acc_ls <- week_accrue(accrual_obj, struct_obj)
      accrual_obj <- week_acc_ls[[1]]
      week_acc <- week_acc_ls[[2]]
      
      # Assign the week's accrual to the object
      if (
        accrual_obj@week <= 
          accrual_obj@target_times[length(accrual_obj@target_times)]
      ) { 
        accrual_obj@accrual[accrual_obj@week, , ] <-
          week_acc
      } else {
        # Extending a predefined array is a pain in the arse in R
        accrual_obj@accrual <- extend_week(accrual_obj@accrual, week_acc)
      }

      # Apply site cap
      accrual_obj <- apply_site_cap(accrual_obj)

      # Apply arm cap, and then adjust site cap appropriately
      obj_list <- apply_arm_cap(accrual_obj, struct_obj)
      accrual_obj <- obj_list[[1]]
      struct_obj <- obj_list[[2]]

    }

    return(list(accrual_obj, struct_obj))
  }


#' Extend predefined array along week axis
#' @param accrual accrual array
#' @param week_acc accrual for week
#' 
extend_week <- function(accrual, week_acc) {
  acc <- aperm(accrual, c(2, 3, 1))
  acc <- array(
    c(
      as.vector(acc),
      as.integer(as.vector(rep(0, length(week_acc))))
    ),
    c(
      dim(acc)[1:2],
      dim(acc)[3] + 1
    )
  ) 

  # Week_acc is arms * sites
  # accrual is weeks * arms * sites permuted to arms * sites * weeks
  acc[, , dim(acc)[3]] <- as.integer(week_acc)
  acc <- aperm(acc, c(3, 1, 2))
  dimnames(acc) <- dimnames(accrual)
  
  acc
}


#' Check whether an object is of class "accrual".
#' 
#' @param x Object to test
#' @param ... For R CMD check compatibility.
#' 
#' @return TRUE or FALSE
#' 
#' @importFrom S7 new_generic method class_any
#' 
#' @export
#' 
is.accrual <- S7::new_generic("is.accrual", "x")
S7::method(is.accrual, S7::class_any) <- function(x) {
  inherits(x, "biomkrAccrual::accrual")
}


