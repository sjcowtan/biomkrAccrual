spine_out <- spine()

test_that("Output of spine() should be an array", {
  checkmate::expect_array(spine_out)
})