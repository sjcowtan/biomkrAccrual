#' A class to contain the accrued patient allocations,
#' by site, treatment/control arm and week
#' 
#' @slot accrual 3-D array with axes site, experimental arm
#' and week
#' @slot accrual_period Number of weeks in planned recruitment period
#' @slot phase_changes Vector of week numbers when arms closed
#' @slot site_closures Vector of weeks sites closed; NA indicates open
#' @slot week Current recruitment week
#' @slot active_arms Vector of indices of open arms
#' @slot active_sites Vector of indices of open sites
#' @slot shared_control TRUE if a shared control arm is being used, else FALSE
#' @slot site_in_region Vector of indices for which set of expected 
#' prevalences each site should use 
#' @slot site_cap Vector of maximum number of patients for each site
#' @slot site_mean_rate Vector of expected recruitment rates for each site
#' @slot site_rate Vector of recruitment rates for each site, drawn from a 
#' gamma distribution
#' @slot start_week Vector of weeks that each site opens recruitment
#' @slot index Vector of index numbers for each site
#' 
#' @param treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @param shared_control TRUE if all experimental arms share one control arm;
#' FALSE if each has their own
#' @param centres_df Dataframe with columns "site", "start_month", "mean_rate", 
#' "region" and "site_cap"
#' @param accrual_period Maximum recruitment period, in weeks
#' 
#' usage accrual(treatment_arm_ids, shared_control, centres_df, accrual_period)
#' @name accrual
#'
#' @export
#' 
#' @importFrom rlang check_dots_empty 
#' 
accrual <- S7::new_class("accrual",
  package = "biomkrAccrual",
  properties = list(
    accrual = S7::class_integer,
    accrual_period = S7::class_integer,
    phase_changes = S7::class_integer,
    site_closures = S7::class_integer,
    week = S7::class_integer,
    active_arms = S7::class_integer,
    active_sites = S7::class_integer,
    shared_control = S7::class_logical,
    site_in_region = S7::class_integer,
    site_cap = S7::class_integer,
    site_mean_rate = S7::class_double,
    site_rate = S7::class_double,
    site_start_week = S7::class_integer,
    site_index = S7::class_integer
  ),
  constructor = function(
    treatment_arm_ids = S7::class_missing,
    shared_control = S7::class_missing,
    centres_df = S7::class_missing,
    accrual_period = S7::class_missing
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
          accrual_period,
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
              paste("C", names(treatment_arm_ids), sep = "-")
            }
          ),
          Centres = c(paste("Centre", unique(centres_df$site)))
        )
      ),
      accrual_period = accrual_period,
      phase_changes = rep(NA_integer_, length(treatment_arm_ids)),
      site_closures = rep(NA_integer_, length(unique(centres_df$site))),
      week = as.integer(1),
      active_arms = seq_len(length(treatment_arm_ids)),
      active_sites = seq_len(length(unique(centres_df$site))),
      shared_control = shared_control,
      site_in_region = as.integer(centres_df$region),
      site_cap = as.integer(centres_df$site_cap),
      site_mean_rate = as.numeric(centres_df$mean_rate),
      site_rate = NA_real_,
      site_start_week = as.integer(centres_df$start_week),
      site_index = as.integer(centres_df$site)
    )
  }
)


#' Sum accrual array by site, to enable checking against site caps
#' @param obj Object of class "accrual"
#' @return vector of total accrual by recruitment site
#' 
site_sums <- S7::new_generic("site_sums", "obj")
S7::method(site_sums, accrual) <- function(obj) {
  # Permute the array so that the first dimension is the
  # dimension you want to get sums for (centre)
  rowSums(aperm(obj@accrual), c(3, 2, 1))
}


#' Sum accrual array by experimental arm (including control)
#' @param obj Object of class "accrual" 
#' @param control_total Logical; if TRUE return single total for control arms
#' @return vector of total accrual by experimental arm
#' 
treat_sums <- S7::new_generic("treat_sums", "obj")
S7::method(treat_sums, accrual) <- function(obj, control_total = FALSE) {
  # Permute the array so that the first dimension is the
  # dimension you want to get sums for (experimental arms)
  arm_sums <-
    rowSums(aperm(obj@accrual, c(2, 1, 3)))

  # If want total for control arms rather than separate values
  if (control_total) {
    # Want individual values for treat arms, but sum of control arms
    if (!obj@shared_control) {
      no_treat_arms <- length(obj@phase_changes)
      arm_sums <- c(
        # Treatment arms
        arm_sums[seq_len(length(obj@phase_changes))],
        # Control
        sum(arm_sums[seq(no_treat_arms + 1, length.out = 2 * no_treat_arms)])
      )
    }
  }
  return(arm_sums)
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
    capped <- sample(population, size = captotal)
  }

  return(capped)
}


#' Get no. weeks from months
#' Currently 4 week month, fix later using lubridate
#' @param months Duration in months
#' @return Duration in weeks
#' 
get_weeks <- function(months) {
  as.integer(round(months * 4, 0))
}


#' Method to increment site rates by gamma-distributed rates of sites
#' opening an accrual pathway in the current week.
#' @param obj An object of type "accrual"
#' @return Modified object with new site rates
#' 
#' @importFrom stats rgamma
#' 
set_site_rates <- S7::new_generic("site_start_rates", "obj")
S7::method(set_site_rates, accrual) <- function(obj, fixed_site_rates) {

  # If this is the first time calling this, initialise with rate 0
  if (any(is.na(obj@site_rate))) {
    obj@site_rate <- rep(0, dim(obj@accrual)[3])
  }

  indices <- which(obj@site_start_week == obj@week)

  if (length(indices) > 0) {

    if (fixed_site_rates) {
      rates <- obj@site_mean_rate(indices) / 4
    } else {
      # mean_rates are in recruitment per month, so scale = 4 converts
      rates <- stats::rgamma(
        n = length(indices),
        shape = obj@site_mean_rate[indices],
        # Per week not per month
        rate = 4
      )
    }

    # When multiple recruitment sources at a given site, want them
    # to stack
    obj@site_rate[obj@site_index[indices]] <- 
      obj@site_rate[obj@site_index[indices]] + rates

  }

  return(obj)
}


#' Implement site cap on a week's accrual. 
#' @param obj Accrual object
#' @param site_cap Maximum number of patients per site
#' @return Modified accrual object with capped week's accrual and 
#' with any capped sites removed from active_sites
#' 
apply_site_cap <- S7::new_generic("apply_site_cap", "obj")
S7::method(apply_site_cap, accrual) <- function(obj) {
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
#' @param target_arm_size Maximum number of patients per arm
#' (can be a vector with a value for each arm, or a scalar)
#' @return Modified accrual object with capped week's accrual
#' 
#' @importFrom rlang abort
#' 
apply_arm_cap <- 
  S7::new_generic("apply_arm_cap", c("accrual_obj", "struct_obj"))
S7::method(apply_arm_cap, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj, target_arm_size) {
    
    # Get totals for experimental arms (dropping control)
    arm_sums <- 
      treat_sums(accrual_obj)[seq_len(length(accrual_obj@phase_changes))]

    # Compare with cap
    arm_captotal <- arm_sums - target_arm_size

    # Inactive arms can be at cap but not exceed it
    if (any(arm_captotal[-accrual_obj@active_arms] > 0)) {
      rlang::abort(paste(
        "Inactive arm exceeded cap:", 
        which(arm_captotal[-accrual_obj@active_arms] > 0)
      ))
    }

    ### Change references to active arms

    # Active arm indices which exceed cap
    active_tocap <- which(arm_captotal > 0)

    if (length(active_tocap) > 0) {
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
    capped_arms <- arm_captotal[accrual_obj@active_arms] >= 0

    accrual_obj@phase_changes[accrual_obj@active_arms[capped_arms]] <- 
      accrual_obj@week
    

    # Update active_arms
    accrual_obj@active_arms <- which(arm_captotal < 0)
    # Also on the trial structure object
    struct_obj <- remove_treat_arms(struct_obj, which(arm_captotal >= 0))

    # Recheck site caps
    site_captotal <- site_sums(accrual_obj) - accrual_obj@site_cap
    
    # Record closing week for capped sites
    capped_sites <- site_captotal[accrual_obj@active_sites] >= 0
    accrual_obj@site_closures[accrual_obj@active_sites[capped_sites]] <-
      accrual_obj@week

    # Update active sites
    accrual_obj@active_sites <- which(site_captotal < 0)

    return(list(accrual_obj, struct_obj)) 
  }




#' Randomise a week's expected accrual amongst the sites, according to 
#' prevalence.
#' @param class_list A list of one object of class "accrual" and one
#' of class "trial_structure".
#' @return Matrix of week's accrual by site and recruitment arm.
#' 
week_accrue <- S7::new_generic("week_accrue", c("accrual_obj", "struct_obj"))
S7::method(week_accrue, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj, fixed_site_rates) {

    # Update the site rates
    accrual_obj <- set_site_rates(accrual_obj, fixed_site_rates)

    # Initialising (sites * experimental arms)
    week_mx <- matrix(
      as.integer(0), 
      nrow = dim(accrual_obj@accrual)[2], 
      ncol = dim(accrual_obj@accrual)[3]
    )

    week_acc <- as.integer(rep(0, dim(accrual_obj@accrual)[3]))

    # Accrual per site is poisson distributed
    week_acc[accrual_obj@active_sites] <- rpois(
      n = length(accrual_obj@active_sites), 
      lambda = accrual_obj@site_rate[accrual_obj@active_sites]
    )

    # Loop over recruiting sites
    for (isite in which(week_acc > 0)) {

      # Total probability for each experimental arm for the 
      # relevant site prevalence set
      probs <- colSums(struct_obj@experimental_arm_prevalence[
        , , accrual_obj@site_in_region[isite]
      ])


      # Sample experimental arms according to probabilities
      assigns <- sample(
        seq_len(length(probs)),
        prob = probs,
        size = week_acc[isite],
        replace = TRUE
      )

      # Increment week's assignment matrix
      for (i in assigns) {
        week_mx[i, isite] <- as.integer(week_mx[i, isite] + 1)
      }
    }

    return(list(accrual_obj, week_mx))
  }


#' Generate accrual for a week and assign it to the 
#' accrual object. Increments the property "week",
#' which holds the next week number to accrue.
#' @param accrual_obj An object of class "accrual"
#' @param struct_obj An object of class "trial_structure"
#' @param target_arm_size Maximum number of patients per arm
#' (can be a vector with a value for each arm, or a scalar)
#' @param fixed_site_rates TRUE if expected site rate to be used; FALSE 
#' draws the site rate from a gamma distribution
#' 
#' @return An object of class "accrual"
#'
accrue_week <- S7::new_generic("accrue_week", c("accrual_obj", "struct_obj"))
S7::method(accrue_week, list(accrual, trial_structure)) <- 
  function(accrual_obj, struct_obj, target_arm_size, fixed_site_rates) {

    # Should not get here if there aren't any but
    if (length(accrual_obj@active_sites) > 0) {

      week_acc_ls <- week_accrue(accrual_obj, struct_obj, fixed_site_rates)
      accrual_obj <- week_acc_ls[[1]]
      week_acc <- week_acc_ls[[2]]

      # Assign the week's accrual to the object
      accrual_obj@accrual[accrual_obj@week, , ] <-
        week_acc#rue(accrual_obj, struct_obj)

      # Apply site cap
      accrual_obj <- apply_site_cap(accrual_obj)

      # Apply arm cap, and then adjust site cap appropriately
      obj_list <- apply_arm_cap(accrual_obj, struct_obj, target_arm_size)
      accrual_obj <- obj_list[[1]]
      struct_obj <- obj_list[[2]]

    }

    return(list(accrual_obj, struct_obj))
  }


