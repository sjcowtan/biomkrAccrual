#' A class for an object to contain the proportions of arriving patients to be
#' allocated to each arm of the trial.
#' 
#' @slot recruit_arm_prevalence Numeric vector of the expected proportion 
#' of patients eligible for each recruitment arm.
#' @slot recruit_arm_names Character vector of the names of the recruitment 
#' arms.
#' @slot treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @slot recruit_arm_id Automatically generated integer vector of the ID 
#' numbers of the recruitment arms.
#' @slot treatment_counts Automatically generated named integer vector of 
#' number of recruitment arms recruiting to each treatment arm.
#' @slot treatment_arm_struct Automatically generated logical matrix of 
#' treatment arms by recruitment arms.
#' @slot treatment_arm_prev Automatically generated matrix of recruitment
#' prevalences of treatment arms by recruitment arms
#' @name trial-structure
#' 
trial_structure <- S7::new_class("trial_structure",
  package = "biomkrAccrual",
  properties = list(
    # These need explicitly setting
    recruit_arm_prevalence = S7::class_numeric,
    recruit_arm_names = S7::class_character,
    treatment_arm_ids = S7::class_list,
    # These are generated from existing properties at the time they execute
    recruit_arm_ids = S7::new_property(
      getter = function(self) seq_len(length(self@recruit_arm_prevalence))
    ),
    treatment_counts = S7::new_property(
      getter = function(self) length(self@treatment_arm_ids)
    ),
    # Will inherit class matrix despite S7 silliness
    treatment_arm_struct = S7::new_property(
      getter = function(self) {
        get_matrix_struct(self@treatment_arm_ids, self@recruit_arm_prevalence)
      }
    ),
    treatment_arm_prev = S7::new_property(
      getter = 
        function(self) {
          get_matrix_prevalence(
            self@treatment_arm_struct, 
            self@recruit_arm_prevalence
          )
        }
    )
  ),
  # Make new instance by calling trial_structure(props_df, arms_ls)
  constructor = function(
    props_df = S7::class_missing, 
    arms_ls = S7::class_missing,
    ...
  ) {
    # Complain if any other arguments given
    rlang::check_dots_empty()

    # Create the object and populate it
    S7::new_object(
      # Parent class
      S7::S7_object(),
      recruit_arm_names = props_df$category, 
      recruit_arm_prevalence = props_df$proportion, 
      treatment_arm_ids = arms_ls
    )
  },
  validator = function(self) {
    if (!is.vector(self@recruit_arm_prevalence)) {
      "Recruitment arm prevalences must be a vector"
    } else if (length(self@recruit_arm_names) != 
        length(self@recruit_arm_prevalence)) {
      "Different number of recruitment arm names supplied than prevalences"
    } #else if (length(self@recruit_arm_prevalence) < 
        #max(unlist(self@treatment_arm_ids))) {
      #"More recruitment arms defined than prevalences specified"
    #}
  }

)

#' Converts trial structure and prevalence information into matrix form
get_matrix_struct <- function(arms_ls, recruit_arm_prevalence) {

  print(recruit_arm_prevalence)

  # Predeclare matrix as no_arms * no_treatments
  no_treats <- length(arms_ls)
  no_recruits <- length(recruit_arm_prevalence)
  
  # This one is logical, to avoid rounding errors
  arm_structure_mx <- 
    matrix(FALSE, max(unlist(
      arms_ls), 
      no_recruits, 
      na.rm = TRUE
    ), no_treats)

  # Loop, changing to TRUE for arms including that treatment
  for (icol in seq_len(no_treats)) {
    
    if (any(!is.na(arms_ls[[icol]]))) {
      arm_structure_mx[unlist(
        arms_ls[[icol]]), 
        icol] <- TRUE
    }
  }

  return(arm_structure_mx)
}

#' Make matrix with the prevalences by treatment arm and recruitment arm
get_matrix_prevalence <- function(arm_structure_mx, recruit_arm_prevalence) {
  arm_prevalence_mx <- 
    matrix(0, nrow = nrow(arm_structure_mx), ncol = ncol(arm_structure_mx))

  # Now loop, replacing 0 with prevalence
  for (irow in seq_len(length(recruit_arm_prevalence))) {
    if (any(arm_structure_mx[irow, ])) {
      arm_prevalence_mx[irow, which(arm_structure_mx[irow, ])] <- 
        recruit_arm_prevalence[irow] / sum(arm_structure_mx[irow, ])
    }
  }
  print(arm_prevalence_mx)

  return(arm_prevalence_mx)
}

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
remove_recruit_arms <- S7::new_generic("remove_recruit_arms", "x")
S7::method(remove_recruit_arms, trial_structure) <- function(x, arms) {

  # Prevalence of removed arms are 0
  x@recruit_arm_prevalence[arms] <- 0
  # Rescale prevalence so sums to 1
  x@recruit_arm_prevalence <- 
    x@recruit_arm_prevalence / sum(x@recruit_arm_prevalence)
  
  # Get existing allocations
  arms_ls <- x@treatment_arm_ids

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
  x@treatment_arm_ids <- arms_ls

  return(x)

}


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
remove_treat_arms <- S7::new_generic("remove_recruit_arms", "x")
S7::method(remove_treat_arms, trial_structure) <- function(x, arms) {
  
  # Remove treatment arms from list
  for (i in arms) {
    x@treatment_arm_ids[[i]] <- NA_integer_
  }


  print(x@treatment_arm_struct)
  
  x@recruit_arm_prevalence[which(colSums(x@treatment_arm_struct) < 1)] <- 0
  x@recruit_arm_prevalence <- 
    x@recruit_arm_prevalence / sum(x@recruit_arm_prevalence)
  print(x@treatment_arm_prev)
  print(x@recruit_arm_prevalence)
  return(x)
}


