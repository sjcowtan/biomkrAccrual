# Testing postscript font name selector

test_that("get_base_family returns a string", {
  checkmate::expect_string(
    get_base_family(),
    null.ok = FALSE,
    na.ok = FALSE
  )
})

test_that("get_base_family returns either Arial or sans", {
  checkmate::expect_choice(
    get_base_family(),
    c("sans", "Arial", "ArialMT")
  )
})