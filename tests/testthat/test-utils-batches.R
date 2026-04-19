# Testing the code for running multiple simulations
output_path <- "../test_biomkrAccrual_output_data/"

biomkrAccrualSim(
  n = 50,
  shared_control = TRUE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  keep_files = TRUE,
  output_path = output_path
)

# And testing it for umbrellas
biomkrAccrualSim(
  n = 50,
  shared_control = FALSE,
  target_times = c(2, 4),
  fixed_centre_starts = TRUE,
  quietly = TRUE,
  keep_files = TRUE,
  output_path = output_path
)

# Clean up by deleting the test data and figure files
unlink(output_path, recursive = TRUE)