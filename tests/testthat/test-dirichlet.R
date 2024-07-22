# Testing rdirichlet_alt
n <- 3
phi <- 2
mu <- c(0.01, 0.1, 0.29, 0.6)

dirichlet_output <- rdirichlet_alt(n, phi, mu)

test_that("Output of rdirichlet_alt is a matrix", {
  expect_s3_class(dirichlet_output, "matrix")
})

testthat("Output of rdirichlet_alt is floating point", {
  expect_type(rdirichlet_alt, "double")
})

test_that("Matrix dimensions from rdirichlet_alt are correct", {
  expect_equal(dim(dirichlet_output), c(n, length(mu)))
})

test_that("Values from rdirichlet_alt are in range (0, 1)", {
  expect_gt(range(x)[1], 0)
  expect_lt(range(x)[2], 1)
})

test_that("Matrix rows from rdirichlet_alt sum to 1", {
  expect_equal(rowSums(dirichlet_output), rep(1, n))
})
