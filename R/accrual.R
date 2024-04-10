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
    active_sites = S7::class_integer,
    shared_control = S7::class_logical
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
      active_sites = seq(sum(centres_df$no_centres)),
      shared_control = shared_control
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
      arm_sums <- c(
        arm_sums[seq_len(length(obj@phase_changes))],
        sum(arm_sums[seq(
          length(obj@phase_changes) + 1, 
          length.out = 2 * length(obj@phase_changes)
        )])
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
