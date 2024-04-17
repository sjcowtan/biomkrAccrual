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
#' @slot treatment_arm_prevalence Automatically generated matrix of recruitment
#' prevalences of treatment arms by recruitment arms
#' @name trial_structure
#' 
trial_structure <- S7::new_class("trial_structure",
  package = "biomkrAccrual",
  properties = list(
    # These need explicitly setting
    recruit_arm_prevalence = S7::class_data.frame,
    recruit_arm_names = S7::class_character,
    treatment_arm_ids = S7::class_list,
    # These are generated from existing properties at the time they execute
    recruit_arm_ids = S7::new_property(
      getter = function(self) seq_len(nrow(self@recruit_arm_prevalence))
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
    # array[recruit arms, treat arms, prevalence set]
    treatment_arm_prevalence = S7::new_property(
      getter = 
        function(self) {
          get_array_prevalence(
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
    centres_df = S7::class_missing,
    ...
  ) {
    # Complain if any other arguments given
    rlang::check_dots_empty()

    # Create the object and populate it
    S7::new_object(
      # Parent class
      S7::S7_object(),
      recruit_arm_names = props_df$category, 
      recruit_arm_prevalence = 
        props_df[, grepl("^proportion_", names(props_df))], 
      treatment_arm_ids = arms_ls
    )
  },
  validator = function(self) {
    if (!is.numeric(self@recruit_arm_prevalence[, -1])) {
      "Recruitment arm prevalences must be numbers"
    } else if (
      length(self@recruit_arm_names) != nrow(self@recruit_arm_prevalence)
    ) {
      "Different number of recruitment arm names supplied than prevalences"
    } else if (
      !(all(sapply(
        self@treatment_arm_ids, 
        function(v) is.integer(unlist(v)) || is.na(v)
      )))
    ) {
      "Elements from the treatment arm list should be integer vectors or NA"
    } else if (
      nrow(self@recruit_arm_prevalence) < 
        max(unlist(self@treatment_arm_ids), na.rm = TRUE)
    ) {
      "More recruitment arms defined than prevalences specified"
    }
  }

)

#' Converts trial structure and prevalence information into matrix form
get_matrix_struct <- function(arms_ls, recruit_arm_prevalence) {

  # Predeclare matrix as no_arms * no_treatments
  no_treats <- length(arms_ls)
  no_recruits <- nrow(recruit_arm_prevalence)
  
  # This one is logical, to avoid rounding errors
  arm_structure_mx <- 
    matrix(FALSE, max(unlist(arms_ls), no_recruits, na.rm = TRUE), no_treats)

  # Loop, changing to TRUE for arms including that treatment
  for (icol in seq_len(no_treats)) {
    
    if (any(!is.na(arms_ls[[icol]]))) {
      arm_structure_mx[unlist(arms_ls[[icol]]), icol] <- TRUE
    }
  }

  return(arm_structure_mx)
}

#' Make array with the prevalences by treatment arm, recruitment arm and 
#' prevalence set
#' @param arm_structure_mx Logical matrix of recruitment arms by treatment arm
#' @param recruit_arm_prevalence Data frame with sets of prevalences, in 
#' columns, for each arm (rows) 
#' @return arm_prevalence_ar Array of prevalences, recruitment arms * 
#' treatment arms * prevalence sets
#' 
get_array_prevalence <- function(arm_structure_mx, recruit_arm_prevalence) {
  arm_prevalence_ar <- 
    array(0, c( 
      nrow(arm_structure_mx), 
      ncol(arm_structure_mx),
      length(recruit_arm_prevalence)
    ))

  # Now loop, replacing 0 with prevalence
  for (iset in seq_len(ncol(recruit_arm_prevalence))) {
    for (irow in seq_len(nrow(recruit_arm_prevalence))) {
      if (any(arm_structure_mx[irow, ])) {
        arm_prevalence_ar[irow, which(arm_structure_mx[irow, ]), iset] <- 
          recruit_arm_prevalence[irow, iset] / sum(arm_structure_mx[irow, ])
      }
    }
  }  

  return(arm_prevalence_ar)
}


#' Close off treatment arms.
#' In the list of treatment arm IDs, replaces the vector of recruitment
#' arm IDs with NA.
#' If any recruitment arms are closed, sets their prevalence to zero and 
#' recalculates prevalence vector, on the assumption that no patients with
#' those characteristics will be recruited from that point.
#' Class object will automatically generate new trial structure and 
#' prevalence matrices.
#' 
#' @param arms vector or scalar of integer arm ID numbers
#' @rdname trial-structure
#' @export
#' 
remove_treat_arms <- S7::new_generic("remove_treat_arms", "obj")
S7::method(remove_treat_arms, trial_structure) <- function(obj, arms) {
  
  # Mark treatment arms as removed using NA; 
  # automatic getter for treatment_arm_struct does the rest
  arms_ls <- obj@treatment_arm_ids

  obj@treatment_arm_ids <- lapply(arms_ls, function(l) {
    # Remove allocations from the recruitment arms in arms
    l <- l[!(l %in% arms)]
    # Remove treatment arm if no more arms recruiting to it
    if (length(l) < 1) {
      return(NA_integer_)
    } else {
      # Return remaining allocations to treatment arm
      return(l)
    }
  })

  # Set prevalence to 0 for any recruitment arms which now have no
  # experimental arms to recruit to
  obj@recruit_arm_prevalence[which(colSums(obj@treatment_arm_struct) < 1)] <- 0
  # Rescale prevalence so it adds to 1
  # (corresponds to recruitment closing for those characteristics)
  obj@recruit_arm_prevalence <- 
    obj@recruit_arm_prevalence / sum(obj@recruit_arm_prevalence)

  return(obj)
}


