###
# !Do not create an R file which comes before this in the alphabet!
###

#' Ensure methods are registered when methods for generics are implemented in 
#' other packages (S3, S4 & S7). This works like export directives in the 
#' NAMESPACE does for S3 and S4.
#' 
#' @name register-S7-methods
.onLoad <- function(...) {
  S7::methods_register()
}

#' Enable usage of <S7_object>@name in package code even for 
#' older versions of R.
#' 
#' @name enable-for-old-R
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL


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
    recruit_arm_id = S7::new_property(
      getter = function(self) seq_len(length(self@recruit_arm_prevalence))
    ),
    treatment_counts = S7::new_property(
      getter = function(self) length(self@treatment_arm_ids)
    ),
    # Will inherit class matrix despite S7 silliness
    treatment_arm_struct = S7::new_property(
      getter = function(self) get_matrix_struct(self@treatment_arm_ids)
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
    } else if (length(self@recruit_arm_prevalence) < 
        max(unlist(self@treatment_arm_ids))) {
      "More recruitment arms defined than prevalences specified"
    }
  }

)
