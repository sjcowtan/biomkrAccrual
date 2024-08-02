### Testing rdirichlet_alt

n <- 3
phi <- 2
mu <- c(0.01, 0.1, 0.29, 0.6)

# Incomplete testing of assertions on input

test_that("Error if n is not integerish", {
  expect_error(rdirichlet_alt(0.5, mu, phi))
})

test_that("Error if mu is not a vector", {
  expect_error(rdirichlet_alt(n, as.matrix(mu, nrow = 4), phi))
})

test_that("Error if mu is of length < 2", {
  expect_error(rdirichlet_alt(n, c(1), phi))
})

test_that("Error if mu is negative", {
  expect_error(rdirichlet_alt(n, c(-0.02, 0.9), phi))
})

test_that("Error if mu is zero", {
  expect_error(rdirichlet_alt(n, c(0.02, 0), phi))
})

test_that("Error if phi is zero", {
  expect_error(rdirichlet_alt(n, mu, 0))
})

test_that("Error if n is not a numeric scalar greater than 0", {
  expect_error(rdirichlet_alt(n, mu, 0))
})


# Testing of successful function execution

dirichlet_output <- rdirichlet_alt(n, mu, phi)

test_that("Output of rdirichlet_alt is a matrix", {
  expect_true(is.matrix(dirichlet_output))
})

test_that("Output of rdirichlet_alt is floating point", {
  expect_type(dirichlet_output, "double")
})

test_that("Matrix dimensions from rdirichlet_alt are correct", {
  expect_equal(dim(dirichlet_output), c(n, length(mu)))
})

test_that("No missing values in output from rdirichlet_alt", {
  expect_equal(sum(is.na(dirichlet_output)), 0)
})

test_that("Values from rdirichlet_alt are greater than 0", {
  expect_gt(range(dirichlet_output)[1], 0)
})

test_that("Values from rdirichlet_alt are less than 1", {
  expect_lt(range(dirichlet_output)[2], 1)
})

test_that("Matrix rows from rdirichlet_alt sum to 1", {
  expect_equal(rowSums(dirichlet_output), rep(1, n))
})


### Testing bio_prevalence
sites_in_region <- sample(1:3, 10, replace = TRUE)
props <- data.frame(
  proportion_1 = c(0.2, 0.5),
  proportion_2 = c(0.35, 0.9),
  proportion_3 = c(0.01, 0.6)
)

test_that("At least one site", {
  expect_error(bio_prevalence(0, props, 10))
})

test_that("Sites must be in regions we have prevalences for", {
  expect_error(bio_prevalence(
    c(sample(1:3, 9, replace = TRUE), nrow(props) + 1), props, 10
  ))
})

test_that("Proportions must not be missing", {
  expect_error(bio_prevalence(
    sites_in_region,
    data.frame(
      proportion_1 = c(0.2, 0.5),
      proportion_2 = c(0.35, NA),
      proportion_3 = c(0.01, 0.6)
    ),
    10
  ))
})

test_that("Proportions must be valid", {
  expect_error(bio_prevalence(
    sites_in_region,
    data.frame(
      proportion_1 = c(0.2, 0.5),
      proportion_2 = c(0.35, 1.6),
      proportion_3 = c(0.01, 0.6)
    ),
    10
  ))
})

bio_prevalence_out <- bio_prevalence(sites_in_region, props, 10)

test_that("bio_prevalence returns a matrix", {
  checkmate::expect_matrix(
    bio_prevalence_out,
    any.missing = FALSE,
    nrows = nrow(props),
    ncols = length(sites_in_region),
    null.ok = FALSE
  )
})

test_that("bio_prevalence matrix contains valid probabilities", {
  checkmate::expect_numeric(
    bio_prevalence_out,
    lower = 0,
    upper = 1
  )
})

test_that("Prevalence sets sum to 1", {
  expect_equal(colSums(bio_prevalence_out), rep(1, length(sites_in_region)))
})
