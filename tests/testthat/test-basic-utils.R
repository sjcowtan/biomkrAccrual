# Testing makeifnot_dir

testpath <- "../biomkrAccrual_output_data"
if (dir.exists("../biomkrAccrual_output_data")) {
  test_path <- paste0(testpath, "/makeadir")
}

test_that("Non-existant directory can be created", {
  expect_no_failure(
    makeifnot_dir(testpath)
  )
  checkmate::expect_directory_exists(
    testpath,
    access = "rwx"
  )
})

