#' A class to contain the accrued patient allocations,
#' by site, treatment/control arm and week
#' 
#' @slot accrual 3-D array with axes site, experimental arm
#' and week
#' @slot phase_changes Vector of week numbers when arms closed
#' @name accrual
#' 
#' @export 
#' 
accrual <- S7::new_class("accrual",
  package = "biomkrAccrual",
  properties = list(
    accrual = S7::class_integer,
    phase_changes = S7::class_integer,
    site_closures = S7::class_integer,
    week = S7::class_integer,
    active_arms = S7::class_integer,
    active_sites = S7::class_integer
  ),
  constructor = function(
    treatment_arm_ids = S7::class_missing,
    shared_control = S7::class_missing,
    centres_df = S7::class_missing,
    accrual_period = S7::class_missing,
    ...
  ) {
    # Complain if any other arguments given
    rlang::check_dots_empty()

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
          sum(centres_df$no_centres)
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
          Centres = c(paste("Centre", seq(sum(centres_df$no_centres))))
        )
      ),
      phase_changes = rep(NA_integer_, length(treatment_arm_ids)),
      site_closures = rep(NA_integer_, length(sum(centres_df$no_centres))),
      week = as.integer(1),
      active_arms = seq_len(length(treatment_arm_ids)),
      active_sites = seq(sum(centres_df$no_centres))
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

  # Total the control arms if required
  if (control_total) {
    # Which arms are control arms?
    control_arms <- 
      grep("^C(-|ontrol)", dimnames(obj@accrual)[[2]], value = FALSE)

    # Want individual values for treat arms, but sum of control arms
    if (length(control_arms) > 1) {
      arm_sums <- c(
        arm_sums[-control_arms],
        sum(arm_sums[control_arms])
      )
    }
  }

  return(arm_sums)
}


#' Select arms/sites to cap
#' @param population
#' @param captotal
#' @return vector of arms to cap (can include multiples of the same arm)
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

#' Implement site cap on a week's accrual. 
#' @param obj Accrual object
#' @param site_cap Maximum number of patients per site
#' @return Modified accrual object with capped week's accrual and 
#' with any capped sites removed from active_sites
#' 
apply_site_cap <- S7::new_generic("apply_site_cap", "obj")
S7::method(apply_site_cap, accrual) <- function(obj, site_cap) {
  # If any sites exceed site cap, remove accrual from randomly
  # selected arms until sites are at cap 
  # Represents sites closing during the week
  site_captotal <- site_sums(obj) - site_cap

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
#' @param obj Accrual object
#' @param target_arm_size Maximum number of patients per site
#' @return Modified accrual object with capped week's accrual
#' 
apply_arm_cap <- S7::new_generic("apply_arm_cap", "obj")
S7::method(apply_arm_cap, accrual) <- function(obj, target_arm_size, site_cap) {

  # Get totals for experimental arms (dropping control)
  arm_sums <- treat_sums(obj)[seq_len(length(obj@phase_changes))]

  # Compare with cap
  arm_captotal <- arm_sums - target_arm_size

  # Active arm indices which exceed cap
  active_tocap <- obj@active_arms[arm_captotal[obj@active_arms] > 0]

  if (length(active_tocap) > 0) {
    for (arm in active_tocap) {
      # Which sites recruited to that arm this week?
      population <- as.integer(unlist(sapply(
        which(obj@accrual[obj@week, arm, ] > 0), 
        function(j) rep(j, obj@accrual[obj@week, arm, j])
      )))

      # Possible accruals to remove
      if (length(population) > 0) {
        capped <- do_choose_cap(population, arm_captotal[arm])
        print(c("Capping arm", arm, "using sites", capped))

        # Remove capped instances
        for (i in capped) {
          obj@accrual[obj@week, arm, i] <- 
            as.integer(obj@accrual[obj@week, arm, i] - 1)
        }
      }
    }
  }

  # Record closing week for capped arms
  obj@phase_changes[obj@active_arms[arm_captotal[obj@active_arms] >= 0]] <- 
    obj@week

  # Update active_arms
  obj@active_arms <- which(arm_captotal < 0)

  # Recheck site caps
  site_captotal <- site_sums(obj) - site_cap

  # Record closing week for capped sites
  obj@site_closures[obj@active_sites[site_captotal[obj@active_sites] >= 0]]

  # Update active sites
  obj@active_sites <- which(site_captotal < 0)

  return(obj)
}


#' Method to generate a week's accrual using an object of class accrual
#' @param obj Instance of class "accrual"
#' @return site * arm matrix of accrual totals for the week
#' 
week_accrue <- S7::new_generic("week_accrue", "obj") 
S7::method(week_accrue, accrual) <- function(obj) {

  # Initialise matrix with no accrual
  week_accrual <- matrix(
    as.integer(0), 
    dim(obj@accrual)[2], 
    dim(obj@accrual)[3]
  )

  # Active arms plus control arms
  arms_to_accrue <- c(
    # Active arms
    obj@active_arms,
    # Control arms 
    seq(length(obj@phase_changes) + 1, dim(obj@accrual)[2])
  )

  # Generate accrual for active sites and active/control arms
  week_accrual[arms_to_accrue, obj@active_sites] <- matrix(
    as.integer(sample(
      1:3, 
      length(obj@active_sites) * length(arms_to_accrue), 
      replace = TRUE
    )), 
    nrow = length(arms_to_accrue), 
    ncol = length(obj@active_sites)
  )

  return(week_accrual)
}


#' Generate accrual for a week and assign it to the 
#' accrual object. Increments the property "week",
#' which holds the next week number to accrue.
#' @param obj An object of class "accrual"
#' 
accrue_week <- S7::new_generic("accrue_week", "obj")
S7::method(accrue_week, accrual) <- function(obj, site_cap, target_arm_size) {

  # For now, just populate randomly
  if (length(obj@active_sites) > 0) {

    # Assign the week's accrual to the object
    obj@accrual[obj@week, , ] <- week_accrue(obj)

    # Apply site cap
    obj <- apply_site_cap(obj, site_cap)

    # Apply arm cap
    obj <- apply_arm_cap(obj, target_arm_size, site_cap)

    # Increment pointer for the next week to accrue
    obj@week <- obj@week + as.integer(1)
  }

  return(obj)
}