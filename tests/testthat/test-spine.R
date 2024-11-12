bma_out <- biomkrAccrual(quietly = TRUE)

test_that("Output of biomkrAccrual() should be an S7 object", {
  checkmate::expect_class(bma_out, "S7_object")
})

test_that("Property accrual of output of biomkrAccrual() should be an array", {
  checkmate::expect_class(bma_out@accrual, "array")
})