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


# Testing get_arm_closures

# Testing accrual_plot_from_file

# Testing accrual_to_long

acc_df <- data.frame(
  T1 = rep(3, 12),
  T2 = rep(5, 12),
  C = rep(8, 12)
)

test_that("accrual_to_long: returns dataframe of correct format", {
  checkmate::expect_data_frame(
    accrual_to_long(acc_df),
    types = c("integerish", "factor"),
    any.missing = FALSE,
    null.ok = FALSE,
    ncols = 3,
    nrows = 36,
    col.names = "named"
  )
})

atl_out <- accrual_to_long(acc_df)

test_that("accrual_to_long: output is of class accrualplotdata", {
  checkmate::expect_class(
    atl_out,
    c("accrualplotdata", "data.frame")
  )
})

test_that("accrual_to_long: correct column names", {
  checkmate::expect_set_equal(
    names(atl_out),
    c("Week", "Arm", "Recruitment"),
    ordered = TRUE
  )
})

test_that("accrual_to_long: Arm is a factor", {
  checkmate::expect_factor(
    atl_out$Arm,
    levels = c("C", "T1", "T2"),
    empty.levels.ok = FALSE,
    null.ok = FALSE
  )
})

test_that("accrual_to_long: Week and Recruitment are integerish", {
  checkmate::expect_integerish(
    atl_out$Week
  )
  checkmate::expect_integerish(
    atl_out$Recruitment
  )
})

test_that("accrual_to_long: Recruitment is cumulatively summed", {
  expect_equal(
    atl_out$Recruitment,
    c(
      seq.int(from = 3, length.out = 12, by = 3),
      seq.int(from = 5, length.out = 12, by = 5),
      seq.int(from = 8, length.out = 12, by = 8)
    )
  )
})

test_that("accrual_to_long: Weeks are consecutive starting at 1", {
  expect_equal(
    atl_out$Week,
    rep(1:12, 3)
  )
})

test_that("accrual_to_long: arm names in expected distribution", {
  expect_equal(
    as.character(atl_out$Arm),
    rep(c("T1", "T2", "C"), each = 12)
  )
})


# Testing plot.accrualplotdata

### Conveniently atl_out is of class accrualplotdata

#### Thanks to Comevussor and hplieninger 
#### https://stackoverflow.com/questions/31038709/
#### how-to-write-a-test-for-a-ggplot-plot

test_that("plot.accrualplotdata: produces an object of class ggplot", {
  expect_silent(
    plot.accrualplotdata(
      atl_out,
      target_arm_size = 40,
      target_control = 80,
      target_interim = 20,
      accrual_period = 12,
      interim_period = 6
    )
  )
  expect_s3_class(
    plot.accrualplotdata(
      atl_out,
      target_arm_size = 40,
      target_control = 80,
      target_interim = 20,
      accrual_period = 12,
      interim_period = 6
    ),
    "ggplot"
  )
})

test_that("plot.accrualplotdata: dispatch is working for this class", {
  expect_silent(
    plot(
      atl_out,
      target_arm_size = 40,
      target_control = 80,
      target_interim = 20,
      accrual_period = 12,
      interim_period = 6
    )
  )
  expect_s3_class(
    plot(
      atl_out,
      target_arm_size = 40,
      target_control = 80,
      target_interim = 20,
      accrual_period = 12,
      interim_period = 6
    ),
    "ggplot"
  )
})

p <- plot.accrualplotdata(
  atl_out,
  target_arm_size = 40,
  target_control = 80,
  target_interim = 20,
  accrual_period = 12,
  interim_period = 6
)

p1 <- plot(
  atl_out,
  target_arm_size = 40,
  target_control = 80,
  target_interim = 20,
  accrual_period = 12,
  interim_period = 6
)
print(names(p))

test_that("plot.accrualplotdata: S3 dispatch works", {
  expect_identical(
    p$coordinates, 
    p1$coordinates
  )
})


test_that("plot.accrualplotdata: Axis labels are correct", {
  expect_identical(p$labels$x, "Week")
  expect_identical(p$labels$y, "Recruitment")
  expect_identical(p1$labels$x, "Week")
  expect_identical(p1$labels$y, "Recruitment")
})

test_that("plot.accrualplotdata: Legend labels are correct", {
  expect_identical(p$labels$group, "Arm")
  expect_identical(p$labels$colour, "Arm")
  expect_identical(p1$labels$group, "Arm")
  expect_identical(p1$labels$colour, "Arm")
})

# Retrieve underlying lists
class(p) <- "list"
class(p1) <- "list"

# Remove unpredictable "environment" element
p$plot_env <- NULL
p1$plot_env <- NULL


