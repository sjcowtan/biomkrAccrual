output_path <- "../test_biomkrAccrual_output_data/"

bma_out <- biomkrAccrual(
  quietly = TRUE, 
  output_path = output_path
)

test_that("Output of biomkrAccrual() should be an S7 object", {
  checkmate::expect_class(bma_out, "S7_object")
})

test_that("Property accrual of output of biomkrAccrual() should be an array", {
  checkmate::expect_class(bma_out@accrual, "array")
})

test_that("Output directory exists or has been created", {
  checkmate::expect_directory_exists(
    output_path,
    access = "rwx"
  )
})

test_that("Figures directory exists or has been created", {
  checkmate::expect_directory_exists(
    paste0(output_path, "figures/"),
    access = "rwx"
  )
})

# Test that data and figure files exist with today's date in test data area
date <- format(Sys.time(), "%Y-%m-%d")

test_that("spine: writes accrual CSV file", {
  checkmate::expect_file_exists(
    file.path(
      substr(output_path, 1, nchar(output_path) - 1),
      list.files(
        path = substr(output_path, 1, nchar(output_path) - 1), 
        pattern = paste0("^accrual-", date, "-.*\\.csv$")
      )
    )
  )
})


test_that("spine: writes at least one PNG file", {
  checkmate::expect_file_exists(
    file.path(
      paste0(output_path, "figures"),
      list.files(
        path = paste0(output_path, "figures"), 
        pattern = paste0("^accrual-", date, "-.*\\.png$")
      )
    )
  )
})


# Clean up by deleting the test data and figure files
unlink(output_path, recursive = TRUE)