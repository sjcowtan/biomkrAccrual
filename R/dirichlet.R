#' Generate a matrix of `n` rows of sets of probabilities 
#' generated from the Dirichlet distribution, each row 
#' summing to 1. Uses the alternative $\mu$, $\phi$ 
#' parameterisation for the Dirichlet distribution, 
#' representing means and precision.
#' 
#' @param n Number of sets of probabilities (defaults to 1)
#' @param mu Vector of mean values for each probability in the set
#' (defaults to c(1, 1, 1)). Must be greater than 0 and finite, and
#' contain at least one value.
#' @param phi Parameter representing precision, where precision is 
#' 1/variance. Must be positive and finite. Defaults to 1.
#'
#' @examples 
#' rdirichlet_alt(n = 3, mu = c(0.001, 0.029, 0.7), phi = 4)
#' 
#' @import checkmate
#' 
rdirichlet_alt <- function(n = 1, mu = c(1, 1, 1), phi = 1) {

  ### Check and convert inputs

  # True if n is "close to an integer"
  checkmate::assert_integerish(
    n, lower = 1, , upper = 10^7, len = 1, any.missing = FALSE
  )

  # mu should be a numeric vector in the range (0, 1)
  checkmate::assert_vector(
    mu, min.len = 1, strict = TRUE, any.missing = FALSE
  )
  checkmate::assert_numeric(
    mu, lower = 0, finite = TRUE
  )

  # phi should be a positive number
  checkmate::assert_number(phi, lower = 10^-7, finite = TRUE)

  # Make n actually an integer
  n <- as.integer(n)

  # Shape parameter (alpha) from mu and phi
  alpha <- phi * c(1, mu)
  no_probs <- length(alpha)

  # Dirichlet is a set of normalised independent gamma(alpha, 1)
  draws <- matrix(
    rgamma(n * no_probs, alpha), 
    ncol = no_probs, 
    byrow = TRUE
  )
  draws <- draws / rowSums(draws)

  return(draws)
}