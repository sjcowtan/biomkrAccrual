# Testing postscript font name selector

test_that("get_base_family returns a string", {
  checkmate::expect_string(
    get_base_family(),
    null.ok = FALSE,
    na.ok = FALSE
  )
})

test_that("get_base_family returns either Arial or sans", {
  expect_match(
    get_base_family(),
    "([Aa]rial)|(sans)",
    perl = TRUE
  )
})