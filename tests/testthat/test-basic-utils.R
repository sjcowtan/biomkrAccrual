# Testing makeifnot_dir

# Put test directory in the test code output folder
testpath <- "../biomkrAccrual_output_data"
if (!dir.exists(testpath)) {
  dir.create(testpath, mode = "0777")
}

testpath <- paste0(testpath, "/makeadir")

test_that("Non-existant directory can be created", {
  expect_success(
    makeifnot_dir(testpath)
  )
  checkmate::expect_directory_exists(
    testpath,
    access = "rwx"
  )
})

test_that("Exits correctly if directory already exists", {
  expect_success(
    makeifnot_dir(testpath)
  )
})

# Make a non-writeable directory and test can't create directory in it
dir.create(paste0(testpath, "/nonwriteable"), mode = "0555")
print("Umask of nonwriteable")
print(Sys.umask(paste0(testpath, "/nonwriteable")))

bad_testpath <- paste0(testpath, "/nonwriteable/bad")

test_that("Fails if parent directory is not writeable", {
  expect_error(
    makeifnot_dir(bad_testpath)
  )
})

# Clean up
unlink(testpath, recursive = TRUE)