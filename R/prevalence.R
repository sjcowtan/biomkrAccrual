#' A class for an object to contain the proportions of arriving patients to be
#' allocated to each arm of the trial.
#' 
#' @slot recruit_arm_prevalence Numeric vector of the expected proportion 
#' of patients eligible for each recruitment arm.
#' @slot recruit_arm_names Character vector of the names of the recruitment 
#' arms.
#' @slot shared_control TRUE if using a shared control arm for all 
#' experimental arms
#' @slot treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @slot recruit_arm_id Automatically generated integer vector of the ID 
#' numbers of the recruitment arms.
#' @slot treatment_counts Automatically generated named integer vector of 
#' number of recruitment arms recruiting to each treatment arm.
#' @slot treatment_arm_struct Automatically generated logical matrix of 
#' treatment arms by recruitment arms.
#' @slot experimental_arm_prevalence Automatically generated matrix of recruitment
#' prevalences of treatment arms by recruitment arms
#' @name trial_structure
#' 
#' @importFrom rlang check_dots_empty
#' 
trial_structure <- S7::new_class("trial_structure",
  package = "biomkrAccrual",
  properties = list(
    # These need explicitly setting
    recruit_arm_prevalence = S7::class_data.frame,
    recruit_arm_names = S7::class_character,
    shared_control = S7::class_logical,
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
    experimental_arm_prevalence = S7::new_property(
      getter = 
        function(self) {
          get_array_prevalence(
            self@treatment_arm_struct, 
            self@recruit_arm_prevalence,
            self@shared_control
          )
        }
    )
  ),
  # Make new instance by calling trial_structure(props_df, arms_ls)
  constructor = function(
    props_df = S7::class_missing, 
    arms_ls = S7::class_missing,
    centres_df = S7::class_missing,
    shared_control = S7::class_missing,
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
      shared_control = shared_control,
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
#' (treatment arms + control arms) * prevalence sets
#' 
get_array_prevalence <- function(arm_structure_mx, recruit_arm_prevalence, shared_control) {
  no_treatments <- ncol(arm_structure_mx)
  no_recruit_arms <- nrow(arm_structure_mx)
  
  arm_prevalence_ar <- 
    array(0, c( 
      nrow(arm_structure_mx), 
      ifelse(
        shared_control, 
        no_treatments + 1,
        no_treatments * 2
      ),
      ncol(recruit_arm_prevalence)
    ))

  # Now loop, replacing 0 with prevalence
  for (iset in seq_len(ncol(recruit_arm_prevalence))) {
    for (irow in seq_len(no_recruit_arms)) {
      # ASSUMING 1:1, half the recruitment goes to the experimental arm
      if (any(arm_structure_mx[irow, ])) {
        # treatment arms
        arm_prevalence_ar[irow, which(arm_structure_mx[irow, ]), iset] <- 
          recruit_arm_prevalence[irow, iset] / 
          (2 * sum(arm_structure_mx[irow, ]))
        # control arms
        if (shared_control) {
          arm_prevalence_ar[irow, no_treatments + 1, iset] <-
            sum(arm_prevalence_ar[irow, seq(no_treatments), iset])
        } else {
          arm_prevalence_ar[
            irow, seq.int(no_treatments + 1, length.out = no_treatments), iset
          ] <- 
            arm_prevalence_ar[irow, seq_len(no_treatments), iset]
        }
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
#' @param arms vector or scalar of integer arm ID numbers to remove
#' 
#' @export
#' 
remove_treat_arms <- S7::new_generic("remove_treat_arms", "obj")
S7::method(remove_treat_arms, trial_structure) <- function(obj, arms) {
  
  # Mark treatment arms as removed using NA; 
  # automatic getter for treatment_arm_struct does the rest
  obj@treatment_arm_ids[arms] <- NA_integer_

  # Set prevalence to 0 for any recruitment arms which now have no
  # experimental arms to recruit to
  obj@recruit_arm_prevalence[which(colSums(obj@treatment_arm_struct) < 1), ] <- 
    0
    
  # Rescale prevalence so it adds to 1 <- not doing this, unrealistic
  # (corresponds to recruitment closing for those characteristics)
  #for (iset in seq_len(ncol(obj@recruit_arm_prevalence))) {
  #  obj@recruit_arm_prevalence[, iset] <- 
  #    obj@recruit_arm_prevalence[, iset] / 
  #    sum(obj@recruit_arm_prevalence[, iset])
  #}

  return(obj)
}


