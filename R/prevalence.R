#' Close off treatment arms.
#' Recalculates prevalence on the assumption that patients for that
#' treatment arm are no longer recruited.
#' Removes treatment arm from the vector of active treatments for each
#' recruitment arm.
#' Class object will automatically generate new trial structure and 
#' prevalence matrices.
#' 
#' @param arms vector or scalar of integer arm ID numbers
#' @rdname trial-structure
#' @export
#' 
remove_recruit_arms <- S7::new_generic("remove_recruit_arms", "arms")
S7::method(remove_recruit_arms, trial_structure) <- function(arms) {

  # Having problems with integers
  #arms <- as.integer(arms)

  # Prevalence of removed arms are 0
  self@recruit_arm_prevalence[arms] <- 0
  # Rescale prevalence so sums to 1
  self@recruit_arm_prevalence <- 
    self@recruit_arm_prevalence / sum(self@recruit_arm_prevalence)
  
  # Get existing allocations
  arms_ls <- self@treatment_arm_ids

  # Loop over treatment arms
  arms_ls <- lapply(arms_ls, function(l) {
    # Remove allocations from the recruitment arms in arms
    l <- l[!(l %in% arms)]
    # Remove treatment arm if no more arms recruiting to it
    if (length(l) < 1) {
      return(NULL)
    # Return remaining allocations to treatment arm
    } else {
      return(l)
    }
  })

  # Set new allocation list on structure object
  self@treatment_arm_ids <- arms_ls

  return(self)

}


#' Converts trial structure and prevalence information into matrix form
get_matrix_struct <- function(arms_ls) {
  # Predeclare matrix as no_arms * no_treatments
  no_treats <- length(arms_ls)
  
  # This one is logical, to avoid rounding errors
  arm_structure_mx <- 
    matrix(FALSE, max(unlist(
      arms_ls), 
      no_treats, 
      na.rm = TRUE
    ), no_treats)

  # Loop, changing to TRUE for arms including that treatment
  for (icol in seq_len(no_treats)) {
    arm_structure_mx[unlist(
      arms_ls[[icol]]), 
      icol] <- TRUE
  }

  return(arm_structure_mx)
}

get_matrix_prevalence <- function(arm_structure_mx, recruit_arm_prevalence) {
  # Now one with the prevalences

  arm_prevalence_mx <- 
    matrix(0, nrow = nrow(arm_structure_mx), ncol = ncol(arm_structure_mx))

  # Now loop, replacing 1 with prevalence
  for (irow in seq_len(length(recruit_arm_prevalence))) {
    arm_prevalence_mx[irow, which(arm_structure_mx[irow, ])] <- 
      arm_prevalence[irow] / sum(arm_structure_mx[irow, ])
  }

  return(arm_prevalence_mx)
}

