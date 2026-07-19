# Testing the code for running multiple simulations
output_path <- "../test_biomkrAccrual_output_data/"

closures_df <- biomkrAccrualSim(
  n = 10,
  var_lambda = 0.25,
  precision = 10,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  output_path = output_path
)

test_that(paste(
  "Dimensions of output correct:", 
  "Shared control, variable site rates and prevalence"
), {
  expect_equal(
    dim(closures_df),
    c(10, 6)
  )
})

test_that(paste(
  "Last row of output correct:", 
  "Shared control, variable site rates and prevalence"
), {
  expect_equal(
    closures_df[10, ],
    data.frame(
      T1 = 5, T2 = 18, T3 = 19, T4 = 16, T5 = 16, Control = 19,
      row.names = as.integer(10)
    )
  )
})

# And testing it for umbrellas
closures_df <- biomkrAccrualSim(
  n = 10,
  var_lambda = 0.25,
  precision = 10,
  shared_control = FALSE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  output_path = output_path
)

test_that("Separate controls, variable site rates and prevalence", {
  expect_equal(
    dim(closures_df),
    c(10, 10)
  )
})

# And for fixed site rates
closures_df <- biomkrAccrualSim(
  n = 10,
  precision = 10,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  fixed_site_rates = TRUE,
  quietly = TRUE,
  output_path = output_path
)

test_that("Shared control, fixed site rates, variable prevalence", {
  expect_equal(
    closures_df[10, ],
    data.frame(
      T1 = 6, T2 = 15, T3 = 20, T4 = 14, T5 = 16, Control = 20,
      row.names = as.integer(10)
    )
  )
})

# And for keep_files = FALSE
closures_df <- biomkrAccrualSim(
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

# keep_files = FALSE keeps only the arm_closures file
test_that("Don't keep files", {
  expect_equal(
    length(list.files(
      path = output_path,
      pattern = "^arm_closures"
    )),
    length(list.files(
      path = output_path,
      pattern = "^arm_totals"
    )) + 1
  )
})


# Testing matrix_to_long

long_df <- matrix_to_long(
  matrix(
    1:12, 
    ncol = 6,
    dimnames = list(NULL, c(paste0("T", 1:5), "Control"))
  )
)

test_that("matrix_to_long produces a dataframe", {
  checkmate::expect_data_frame(
    long_df,
    ncols = 3,
    nrows = 12,
    types = c("character", "integerish"),
    any.missing = FALSE
  )
})

test_that("matrix_to_long has correct columns", {
  expect_equal(
    names(long_df),
    c("Arm", "Recruitment", "Run")
  )
})

test_that("matrix_to_long contains correct values", {
  expect_equal(long_df$Recruitment, 1:12)
})

# Clean up by deleting the test data and figure files
unlink(output_path, recursive = TRUE)