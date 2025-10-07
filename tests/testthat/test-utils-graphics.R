# Testing theme creation

test_that("theme_bma creates a valid theme", {
  expect_list(theme_bma())
  checkmate::expect_class(theme_bma(), c("theme", "gg"))
})


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


# Testing 