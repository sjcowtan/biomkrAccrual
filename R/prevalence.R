#' A class for an object to contain the proportions of arriving patients to be
#' allocated to each arm of the trial.
#' 
#' @slot recruit_arm_prevalence Numeric vector of the expected proportion 
#' of patients eligible for each recruitment arm.
#' @slot recruit_arm_prevalence_start Numeric vector of the initial proportions
#' of patients eligible for each recruitment arm.
#' @slot recruit_arm_names Character vector of the names of the recruitment 
#' arms.
#' @slot shared_control TRUE if using a shared control arm for all 
#' experimental arms.
#' @slot control_ratio Proportion of patients assigned to control
#' @slot treatment_arm_ids Named list of lists of recruitment arms by 
#' treatment arm.
#' @slot treatment_arm_ids_start Named list of lists of the initial 
#' configuration of recruitment arms by treatment arm.
#' @slot recruit_arm_id Automatically generated integer vector of the ID 
#' numbers of the recruitment arms.
#' @slot treatment_counts Automatically generated named integer vector of 
#' number of recruitment arms recruiting to each treatment arm.
#' @slot treatment_arm_struct Automatically generated logical matrix of 
#' treatment arms by recruitment arms.
#' @slot treatment_arm_struct_start Automatically generated logical matrix of 
#' the initial configuration of treatment arms by recruitment arms.
#' @slot experimental_arm_prevalence Automatically generated matrix of 
#' prevalences of treatment arms by recruitment arms
#' 
#' @param props_df Dataframe of expected biomarker prevalences for the 
#' regions, with one column `category` containing names for the 
#' biomarkers, and one column per region.
#' @param arms_ls List of lists of recruitment arms which recruit to
#' each treatment arm.
#' @param centres_df Dataframe containing columns `site`, the index
#' number of each site; `start_month`, the month in which that site
#' is expected to start recruiting; `mean_rate`, the expected number 
#' of patients from that site per month; `region`, the index of the
#' region the site is in (should be in the same order as the columns
#' in `props_df`); and an optional column `site_cap`, if there is a 
#' recruitment cap on any of the sites.
#' @param precision For the Dirichlet model of biomarker prevalences, 
#' variability decreases as precision increases. Defaults to 10.
#' @param shared_control TRUE if all experimental arms share one 
#' control arm; FALSE if they each have separate control arms.
#' @param control_ratio Proportion of patients assigned to control
#' @param fixed_region_prevalences TRUE if biomarker prevalences 
#' should be considered to be identical for all sites within a 
#' region; FALSE if they should be drawn from a Dirichlet distribution
#' with a mean of the specified prevalence.
#' 
#' @name trial_structure
#' 
#' @import S7
#' 
trial_structure <- S7::new_class("trial_structure",
  package = "biomkrAccrual",
  properties = list(
    # These need explicitly setting
    recruit_arm_prevalence = S7::class_double,
    recruit_arm_prevalence_start = S7::class_double,
    recruit_arm_names = S7::class_character,
    shared_control = S7::class_logical,
    control_ratio = S7::class_vector,
    treatment_arm_ids = S7::class_list,
    treatment_arm_ids_start = S7::class_list,
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
    treatment_arm_struct_start = S7::new_property(
      getter = function(self) {
        get_matrix_struct(
          self@treatment_arm_ids_start, self@recruit_arm_prevalence_start
        )
      }
    ),
    # array[recruit arms, treat arms, prevalence set]
    experimental_arm_prevalence = S7::new_property(
      getter = 
        function(self) {
          get_array_prevalence(
            self@treatment_arm_struct, 
            self@recruit_arm_prevalence,
            self@shared_control,
            self@control_ratio
          )
        }
    )
  ),
  # Make new instance by calling trial_structure(props_df, arms_ls)
  constructor = function(
    props_df = S7::class_missing, 
    arms_ls = S7::class_missing,
    centres_df = S7::class_missing,
    precision = S7::class_missing,
    shared_control = S7::class_missing,
    control_ratio = S7::class_missing,
    fixed_region_prevalences = S7::class_missing
  ) {
    # Create the object and populate it
    S7::new_object(
      # Parent class
      S7::S7_object(),
      recruit_arm_names = props_df$category, 
      recruit_arm_prevalence = 
        get_recruit_arm_prevalence(
          props_df, centres_df, precision, fixed_region_prevalences
        ),
      recruit_arm_prevalence_start = 
        get_recruit_arm_prevalence(
          props_df, centres_df, precision, fixed_region_prevalences
        ),
      shared_control = shared_control,
      control_ratio = control_ratio,
      treatment_arm_ids = arms_ls,
      treatment_arm_ids_start = arms_ls
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


#' Dirichet regression model for biomarker prevalence
#' 
#' For each site, draws from the Dirichlet distribution using the
#' expected prevalences for the region. 
#' expected prevalences by region
#' 
#' @param props_df Dataframe with one row per biomarker.  Has one column 
#' `category`, containing names of the biomarkers, plus one column for 
#' each region, containing the expected biomarker prevalences for the regions.
#' @param centres_df Dataframe with one row per site, including one column
#' `region`, containing the index number for the region for each site.
#' Indices are assumed to be in the same order as the columns in `props_df`.
#' @param precision Variability decreases as precision increases.
#' @param fixed_region_prevalences TRUE if biomarker prevalences 
#' should be considered to be identical for all sites within a 
#' region; FALSE if they should be drawn from a Dirichlet distribution
#' with a mean of the specified prevalence.
#'   
#' @return Matrix of prevalences with one column per site and one 
#' row per biomarker; each column sums to 1.
#' 
#' @importFrom checkmate assert_names assert_data_frame assert_atomic_vector 
#' assert_integerish assert_numeric
#' 
get_recruit_arm_prevalence <- function(
  props_df, centres_df, precision, fixed_region_prevalences
) {

  # Check format of fixed_region_prevalences
  checkmate::assert_logical(
    fixed_region_prevalences, 
    len = 1, 
    any.missing = FALSE, 
    null.ok = FALSE
  )

  # If fixed_region_prevalences is FALSE, check contents of 
  # precision; otherwise should be NULL
  if (fixed_region_prevalences) {
    assert_null(precision)
  } else {
    # Check format and content of precision
    checkmate::assert_numeric(
      precision,
      lower = 10^-7,
      finite = TRUE,
      len = 1,
      any.missing = FALSE,
      null.ok = FALSE
    )
  }

  # Check format and content of centres_df

  checkmate::assert_data_frame(
    centres_df,
    types = "numeric",
    any.missing = FALSE,
    min.cols = 5,
    max.cols = 6,
    min.rows = 1,
    col.names = "named",
    null.ok = FALSE
  )

  checkmate::assert_names(
    names(centres_df),
    subset.of = c(
      "site", "start_month", "mean_rate", "region", "site_cap", "start_week"
    ),
    must.include = c(
      "site", "start_month", "mean_rate", "region", "start_week"
    )
  )

  sites_in_region <- centres_df$region

  checkmate::assert_atomic_vector(
    sites_in_region,
    any.missing = FALSE,
    min.len = 1,
    max.len = 10^4
  )

  checkmate::assert_integerish(
    sites_in_region,
    lower = 1,
    upper = 10^4,
    any.missing = FALSE,
    null.ok = FALSE
  )
  
  # Check format and content of props_df

  checkmate::assert_data_frame(
    props_df,
    types = c("numeric", "character"),
    any.missing = FALSE,
    # number of regions + "category"
    min.cols = max(sites_in_region) + 1,
    min.rows = 2,
    null.ok = FALSE
  )

  checkmate::assert_names(names(props_df), must.include = "category")


  # Any column that isn't "category" is assumed to be a region -
  # allows for named regions
  region_prevalence <- 
    props_df[, grep("^category$", names(props_df), invert = TRUE)]

  checkmate::assert_numeric(
    as.matrix(region_prevalence),
    lower = 0,
    upper = 1,
    finite = TRUE
  )

  if (fixed_region_prevalences) {
    # Use region prevalences unchanged
    recruit_arm_prevalence <- as.matrix(
      region_prevalence[, sites_in_region]
    )
    # Scale columns to sum to 1
    recruit_arm_prevalence <- sweep(
      recruit_arm_prevalence,
      2,
      colSums(recruit_arm_prevalence),
      FUN = "/"
    )

  } else {
    
    # Draw from Dirichlet distribution
    recruit_arm_prevalence <- do_dirichlet_draws(
      region_prevalence, sites_in_region, precision
    )
  }

  return(recruit_arm_prevalence)

}


#' Draws from dirichlet regression model for biomarker prevalences
#' to create the prevalence matrix.
#' 
#' @param region_prevalence Dataframe with one column for each 
#' region and one row for each biomarker, containing prevalences as 
#' probabilities.
#' @param sites_in_region Vector with a region index number for each
#' site, the index determined by the order of the columns in
#' `region_prevalence`.
#' @param precision Variability decreases as precision increases.
#'  
do_dirichlet_draws <- function(region_prevalence, sites_in_region, precision) {
  
  # Predeclare prevalence matrix
  recruit_arm_prevalence_mx <- matrix(
    data = NA, 
    ncol = length(sites_in_region),
    nrow = nrow(region_prevalence)
  )
  
  # Draw prevalences for sites in each region in turn
  for (region in unique(sites_in_region)) {

    # Sites in region
    site_indices <- which(sites_in_region == region)

    bio_prevalence <- rdirichlet_alt(
      n = length(site_indices),
      mu = region_prevalence[, region],
      phi = precision
    )
    # rdirichlet_alt produces the transpose of what we want,
    # for consistency with other implementations of rdirichlet
    recruit_arm_prevalence_mx[, site_indices] <- t(bio_prevalence)
  }

  return(recruit_arm_prevalence_mx)
}


#' Converts trial structure and prevalence information into matrix form
#' 
#' @param arms_ls List of lists of recruitment arms which recruit to
#' each treatment arm.
#' @param recruit_arm_prevalence Matrix of prevalences with one row per 
#' biomarker and one column per site; each column sums to 1.
#' 
#' @return Logical matrix with one row per biomarker and one column per
#' treatment arm; TRUE where the treatment recruits from that biomarker.
#' 
get_matrix_struct <- function(arms_ls, recruit_arm_prevalence) {

  # Predeclare matrix as no_biomarkers * no_treatments
  no_treats <- length(arms_ls)
  no_biomarkers <- nrow(recruit_arm_prevalence)
  
  # This one is logical, to avoid rounding errors
  arm_structure_mx <- 
    matrix(FALSE, max(unlist(arms_ls), no_biomarkers, na.rm = TRUE), no_treats)

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
#' @param shared_control TRUE if all treatment arms share the
#' same control arm; FALSE if each treatment arm has its own 
#' control. Defaults to TRUE.
#' @param control_ratio Ratio of patients assigned to treatment versus control
#' 
#' @return arm_prevalence_ar Array of prevalences, recruitment arms * 
#' (treatment arms + control arms) * prevalence sets
#' 
get_array_prevalence <- function(
  arm_structure_mx, recruit_arm_prevalence, shared_control, control_ratio
) {
  no_treatments <- ncol(arm_structure_mx)
  no_recruit_arms <- nrow(arm_structure_mx)
  no_regions <- ncol(recruit_arm_prevalence)

  scale_factor <- 
    ifelse(
      sapply(
        rowSums(arm_structure_mx), 
        function(x) isTRUE(all.equal(x, 0, tolerance = 1e-6))
      ),
      1, 
      rowSums(arm_structure_mx)
    )

  # List of matrices by centre 
  ## Total proportion recruited for centre on each open arm 
  ## - more efficient to fix later
  prev_ls <- lapply(
    seq(no_regions),
    function(i) recruit_arm_prevalence[, i] * arm_structure_mx / scale_factor
  )

  # Apply control configuration
  if (shared_control) {
    prev_ls <- lapply(
      prev_ls,
      function(mx) cbind(mx * control_ratio[1], rowSums(mx) * control_ratio[2])
    )
  } else {
    prev_ls <- lapply(
      prev_ls,
      function(mx) cbind(mx * control_ratio[1], mx * control_ratio[2])
    )
  }

  # Make into array
  arm_prevalence_ar <- array(
    data = do.call(cbind, prev_ls),
    dim = c(dim(prev_ls[[1]]), length(prev_ls))
  )
  
  return(arm_prevalence_ar)
}


#' Close off treatment arms.
#' In the list of treatment arm IDs, replaces the vector of recruitment
#' arm IDs with NA.
#' If any recruitment arms are closed, sets their prevalence to zero but 
#' does not recalculate prevalence vector, on the assumption that recruitment
#' from aites will fall because no patients with those characteristics will 
#' be recruited from that point.
#' Class object will automatically generate new trial structure and 
#' prevalence matrices.
#' 
#' @param structure_obj An object of class `trial_structure`.
#' @param arms Vector or scalar of integer arm ID numbers to remove
#' 
#' @usage remove_treat_arms(structure_obj, arms)
#' 
remove_treat_arms <- S7::new_generic("remove_treat_arms", "structure_obj")
S7::method(remove_treat_arms, trial_structure) <- function(
  structure_obj,
  arms
) {
  
  # Mark treatment arms as removed using NA; 
  # automatic getter for treatment_arm_struct does the rest
  structure_obj@treatment_arm_ids[arms] <- NA_integer_

  # Set prevalence to 0 for any recruitment arms which now have no
  # experimental arms to recruit to
  no_exp_arms <- rowSums(structure_obj@treatment_arm_struct) < 1
  structure_obj@recruit_arm_prevalence[no_exp_arms, ] <- 0

  return(structure_obj)
}


#' Check whether an object is of class "trial_structure".
#' 
#' @param x Object to test
#' 
#' @return TRUE or FALSE
#' 
#' @importFrom S7 new_generic method class_any
#' @export
#' 
is.trial_structure <- S7::new_generic("is.trial_structure", "x")
S7::method(is.trial_structure, S7::class_any) <- function(x) {
  inherits(x, "biomkrAccrual::trial_structure")
}

