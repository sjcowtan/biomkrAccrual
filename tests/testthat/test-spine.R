spine_out <- spine()

test_that("Output of spine() should be a list of length 2", {
  checkmate::expect_list(spine_out, len = 2)
})

test_that("Second element of spine_output should be an integer vector", {
  checkmate::expect_atomic_vector(
    spine_out[[2]], names = "named", min.len = 2
  )
})