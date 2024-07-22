#' Generate a matrix of `n` rows of sets of probabilities 
#' generated from the Dirichlet distribution, each row 
#' summing to 1. Uses the alternative $\mu$, $\phi$ 
#' parameterisation for the Dirichlet distribution, 
#' representing means and precision.
#' 
#' @param n Number of sets of probabilities (defaults to 1)
#' @param mu Vector of mean values for each probability in the set
#' @param phi Parameter representing precision, where precision is 
#' 1/variance. Must be positive.
#'
#' @example rdirichlet_alt(3, c(0.001, 0.029, 0.7), 4)
rdirichlet_alt <- function(n = 1, mu = c(0.2, 0.3, 0.5), phi = 1) {

  ### Check and convert inputs

  # True if n is "close to an integer"
  checkmate::assert_integerish(n, lower = 1, len = 1)

  # mu should be a numeric vector in the range (0, 1)
  checkmate::assert_vector(
    mu, min.len = 2, strict = TRUE, any.missing = FALSE
  )
  checkmate::assert_numeric(
    mu, lower = 10^-10, upper = 1 - 10^-10
  )

  # phi should be a positive number
  checkmate::assert_number(phi, lower = 10^-10)

  # Make n actually an integer
  n <- as.integer(n)

}