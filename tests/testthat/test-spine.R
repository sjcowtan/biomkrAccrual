bma_out <- biomkrAccrual(quietly = TRUE)

test_that("Output of biomkrAccrual() should be an S7 object", {
  checkmate::expect_class(bma_out, "S7_object")
})

test_that("Property accrual of output of biomkrAccrual() should be an array", {
  checkmate::expect_class(bma_out@accrual, "array")
})

test_that("Output directory exists or has been created", {
  checkmate::expect_directory_exists(
    "../biomkrAccrual_output_data/",
    access = "rwx"
  )
})

test_that("Figures directory exists or has been created", {
  checkmate::expect_directory_exists(
    "../biomkrAccrual_output_data/figures/",
    access = "rwx"
  )
})

# Test that data and figure files exist with today's date in test data area
date <- format(Sys.time(), "%y-%m-%d")
print(getwd())

## Can't put wildcard in test
files <- list.files(
  "../biomkrAccrual_output_data/",
  pattern = "accrual-", date, "_*\\.csv"
)
print(files)

test_that("spine: writes accrual CSV file", {
  checkmate::expect_file_exists(
    paste0("../biomkrAccrual_output_data/accrual-", date, "_*.csv")
  )
})


# Clean up by deleting the data and figure files with today's date