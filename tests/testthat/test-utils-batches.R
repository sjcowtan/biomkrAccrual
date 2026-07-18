# Testing the code for running multiple simulations
output_path <- "../test_biomkrAccrual_output_data/"

closures_mx <- biomkrAccrualSim(
  n = 50,
  var_lambda = 0.25,
  precision = 10,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  output_path = output_path
)

test_that("Shared control, variable site rates and prevalence", {
  expect_equal(
    dim(closures_mx),
    c(50, 6)
  )
})

# And testing it for umbrellas
closures_mx <- biomkrAccrualSim(
  n = 10,
  var_lambda = 0.25,
  precision = 10,
  shared_control = FALSE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  output_path = output_path
)

print(head(closures_mx))

test_that("Separate controls, variable site rates and prevalence", {
  expect_equal(
    dim(closures_mx),
    c(10, 10)
  )
})

# And for fixed site rates
closures_mx <- biomkrAccrualSim(
  n = 10,
  precision = 10,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  fixed_site_rates = TRUE,
  quietly = TRUE,
  output_path = output_path
)

print((closures_mx))
print(dim(closures_mx))
test_that("Shared control, fixed site rates, variable prevalence", {
  expect_equal(
    dim(closures_mx),
    c(10, 6)
  )
})

# And for keep_files = FALSE
closures_mx <- biomkrAccrualSim(
  n = 10,
  precision = 10,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  fixed_site_rates = TRUE,
  quietly = TRUE,
  keep_files = FALSE,
  output_path = output_path
)

# Clean up by deleting the test data and figure files
unlink(output_path, recursive = TRUE)
print(output_path)

test_that("Don't keep files", {
  expect_equal(
    dim(closures_mx),
    c(10, 6)
  )
})

print(list.files(path = output_path), recursive = TRUE)

# Clean up by deleting the test data and figure files
#unlink(output_path, recursive = TRUE)