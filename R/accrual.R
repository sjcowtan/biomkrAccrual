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
    week = S7::class_integer,
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
      week = as.integer(1),
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
  # dimension you want to get sums for (experimenta arm)
  arm_sums <- rowSums(aperm(obj@accrual), c(2, 1, 3))

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


#' Generate accrual for a week and assign it to the 
#' accrual object. Increments the property "week",
#' which holds the next week number to accrue.
#' @param obj An object of class "accrual"
#' 
accrue_week <- S7::new_generic("accrue_week", "obj")
S7::method(accrue_week, accrual) <- function(obj, site_cap, target_arm_size) {

  # For now, just add 1 to everything
  if (length(obj@active_sites) > 0) {
    week_accrual <- matrix(
      as.integer(sample(
        0:2, 
        dim(obj@accrual)[2] * dim(obj@accrual)[3], 
        replace = TRUE)
      ), 
      nrow = dim(obj@accrual)[2], 
      ncol = dim(obj@accrual)[3]
    )
    week_accrual[, -obj@active_sites] <- as.integer(0)
  }

  print(week_accrual)

  # Assign the week's accrual to the object
  obj@accrual[obj@week, , ] <- week_accrual

  # If any sites exceed site cap, remove accrual from randomly
  # selected arms until sites are at cap 
  # Represents sites closing during the week
  site_captotal <- site_sums(obj) - site_cap
  print(site_captotal)

  if (any(site_captotal > 0)) {
    for (site in which(site_captotal > 0)) {
      # One at a time so don't take too much from one arm
      for (i in seq(site_captotal[site])) {
        # Pick an arm that still accrued this week
        print(obj@accrual[obj@week, , site])
        arm <- sample(which(obj@accrual[obj@week, , site] > 0))
        print(paste("Drop from arm", arm, "total", obj@accrual[obj@week, arm, site]))
        obj@accrual[obj@week, arm, site] <- 
          as.integer(obj@accrual[obj@week, arm, site] - 1)
      }
    }
  }

  # Take full sites out of recruitment
  obj@active_sites <- which(site_captotal >= 0)


  # Increment pointer for the next week to accrue
  obj@week <- obj@week + as.integer(1)

  return(obj)

}