#' Generate a matrix of `n` rows of sets of probabilities 
#' generated from the Dirichlet distribution, each row 
#' summing to 1. Uses the alternative \eqn{\mu} \eqn{\phi} 
#' parameterisation for the Dirichlet distribution, 
#' representing means and precision.
#' 
#' @param n Number of sets of probabilities (defaults to 1)
#' @param mu Vector of mean values for each probability in the set
#' (defaults to c(0.001, 0.029. 0.7)). Must be greater than 0 and 
#' finite, and contain at least two values.
#' @param phi Parameter representing precision, where precision is 
#' 1/variance. Must be positive and finite. Defaults to 10.
#'
#' @example rdirichlet_alt(n = 3, mu = c(0.001, 0.029, 0.7), phi = 10)
#' 
#' @export 
#' 
#' @import checkmate
#' @importFrom stats rgamma
#' 
rdirichlet_alt <- function(
  n = 1, 
  mu = c(0.3, 0.3, 0.3), 
  phi = 10
) {

  ### Check and convert inputs

  # True if n is "close to an integer"
  checkmate::assert_integerish(
    n, lower = 1, , upper = 10^7, len = 1, any.missing = FALSE
  )

  # mu should be a numeric vector in the range [0, 1)
  checkmate::assert_vector(
    mu, min.len = 2, strict = TRUE, any.missing = FALSE
  )
  checkmate::assert_numeric(
    mu, lower = 10^-7, upper = 1 - 10^-7
  )

  # phi should be a positive number
  checkmate::assert_number(phi, lower = 10^-7, finite = TRUE)

  # Make n actually an integer
  n <- as.integer(n)

  # Shape parameter (alpha) from mu and phi; alpha_0 = phi
  alpha <- phi * mu / mu[1]
  no_probs <- length(mu)

  # Dirichlet is a set of normalised independent gamma(alpha, 1)
  draws_mx <- matrix(
    stats::rgamma(n * no_probs, alpha), 
    ncol = no_probs, 
    byrow = TRUE
  )
  draws_mx <- draws_mx / rowSums(draws_mx)

  return(draws_mx)
}


#' Dirichet regression model for biomarker prevalence
#' 
#' For each site, draws from the Dirichlet distribution using the
#' expected prevalences for the region. 
#' 
#' @param site_in_region Integer vector of region IDs (matching the 
#' column numbers in `recruit_arm_prevalence`) for each site
#' @param recruit_arm_prevalence Dataframe of expected prevalences,
#' one column per region, one row per biomarker
#' @param precision Measure of variability (analogous to 1 / variance)
#' 
#' @return Matrix of prevalences with one column per site and one 
#' row per biomarker; each column sums to 1
#' 
bio_prevalence <- function(
  site_in_region,
  recruit_arm_prevalence,
  precision = 10
) {

  # Predeclare prevalence matrix
  bio_prevalences <- matrix(
    data = NA, 
    ncol = length(site_in_region),
    nrow = nrow(recruit_arm_prevalence)
  )

  for (region in unique(site_in_region)) {

    # Sites in region
    positions <- which(site_in_region == region)

    bio_prevalence <- rdirichlet_alt(
      length(positions),
      recruit_arm_prevalence[, region],
      precision
    )

    bio_prevalences[, positions] <- t(bio_prevalence)
  }
  return(bio_prevalences)
}