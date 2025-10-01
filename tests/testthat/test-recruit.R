### Testing accrual() constructor

fixed_acc_obj <- accrual(
  treatment_arm_ids = list(T1 = as.integer(1), T2 = as.integer(2)),
  shared_control = TRUE,
  accrual_period = as.integer(12),
  interim_period = as.integer(6),
  control_ratio = c(1, 1),
  fixed_site_rates = TRUE,
  var_lambda = 0.25,
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 2),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 5)
  )
)

# Testing constructor

test_that(paste(
  "accrual constructor: produces an object of classes",
  "S7_object and biomkrAccrual::accrual"
), {
  checkmate::expect_class(fixed_acc_obj, c(
    "biomkrAccrual::accrual",
    "S7_object"
  ))
})


### Testing is.accrual()

test_that(paste(
  "is.accrual: recognises an object of class",
  "biomkrAccrual::accrual"
), {
  expect_true(
    is.accrual(fixed_acc_obj)
  )
})


### Testing treat_sums()
arr <- array(1:24, 2:4)

test_that("treat_sums: sum by treatment a 3-D accrual array", {
  checkmate::expect_integer(
    treat_sums(arr),
    min.len = 1,
    max.len = ncol(arr),
    lower = 0,
    any.missing = FALSE,
    null.ok = FALSE
  )
})

test_that("treat_sums: output correct for valid array", {
  expect_identical(
    treat_sums(arr),
    as.integer(c(84, 100, 116))
  )
})

## Testing treat_sums on an accrual object

test_that("treat_sums: output correct  for valid accrual object", {
  expect_identical(
    treat_sums(fixed_acc_obj),
    as.integer(c(0, 0, 0))
  )
})

## Testing capping
uncapped <- c(1, 1, 2, 5, 8)
capped <- do_choose_cap(uncapped, 3)

test_that("do_choose_cap: caps population correctly", {
  expect_identical(
    length(capped),
    as.integer(3)
  )
})

# Need to be able to compare frequencies for elements in capped only
table_uncapped <- as.data.frame(table(uncapped))
table_capped <- as.data.frame(table(capped))

test_that("do_choose_cap: does not produce extra repeats", {
  expect_true(
    all(
      table_uncapped[
        which(table_uncapped$uncapped %in% table_capped$capped),
      ]$Freq >= table_capped$Freq
    )
  )
})



# Testing set_site_rates()

set.seed(123)

rates <- set_site_rates(fixed_acc_obj)

test_that("set_site_rates: fixed site rate is correct", {
  expect_equal(rates@site_rate, c(2.5, 0), tolerance = 1e-6)
})


variable_acc_obj <- accrual(
  treatment_arm_ids = list(T1 = as.integer(1), T2 = as.integer(2)),
  shared_control = TRUE,
  accrual_period = as.integer(12),
  interim_period = as.integer(6),
  control_ratio = c(1, 1),
  fixed_site_rates = FALSE,
  var_lambda = 0.25,
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 2),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 5)
  )
)

set.seed(123)

rates <- set_site_rates(variable_acc_obj)

test_that("set_site_rates: variable site rate is correct", {
  expect_equal(rates@site_rate, c(2.427350, 0.0), tolerance = 1e-6)
})


# Testing week_accrue()

## Need a structure object
ts_obj <- trial_structure(
  props_df =  data.frame(
    category = c("B1", "B2", "B3"),
    region_1 = 0.031, 0.454, 0.515
  ),
  arms_ls = list(
    T1 = 1:2,
    T2 = 2:3
  ),
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 2),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 4)
  ),
  precision = 10,
  shared_control = TRUE,
  control_ratio = c(1, 1),
  fixed_region_prevalences = FALSE
)

set.seed(123)

wa_out_ls <- week_accrue(variable_acc_obj, ts_obj)

test_that("week_accrue: first output is an accrual object", {
  expect_equal(
    class(wa_out_ls[[1]]),
    c("biomkrAccrual::accrual", "S7_object")
  )
})

test_that("week_accrue: second output is a valid accrual matrix", {
  checkmate::expect_matrix(
    wa_out_ls[[2]],
    any.missing = FALSE,
    nrows = 3,
    ncols = 2,
    null.ok = FALSE,
    mode = "integer"
  )
})

test_that("week_accrue: accrual values as expected", {
  expect_equal(
    as.vector(wa_out_ls[[2]]),
    c(1, 2, 1, 0, 0, 0)
  )
})

test_that("week_accrue: variable site rate is correct", {
  expect_equal(wa_out_ls[[1]]@site_rate, c(2.427350, 0.0), tolerance = 1e-6)
})


# Testing accrue_week

set.seed(123)

aw_out_ls <- accrue_week(variable_acc_obj, ts_obj)

test_that("accrue_week: first output is an accrual object", {
  expect_equal(
    class(aw_out_ls[[1]]),
    c("biomkrAccrual::accrual", "S7_object")
  )
})

test_that("accrue_week: first output is a structure object", {
  expect_equal(
    class(aw_out_ls[[2]]),
    c("biomkrAccrual::trial_structure", "S7_object")
  )
})

test_that("accrue_week: recruitment totals are a 3D integer array", {
  checkmate::expect_array(
    aw_out_ls[[1]]@accrual,
    mode = "integer",
    any.missing = FALSE,
    d = 3,
    null.ok = FALSE
  )
})

test_that("accrue_week: recruitment totals are correct", {
  expect_equal(
    as.vector(aw_out_ls[[1]]@accrual), 
    c(1, rep(0, 11), 2, rep(0, 11), 1, rep(0, 47))
  )
})

test_that("accrue_week: variable site rate is correct", {
  expect_equal(
    aw_out_ls[[1]]@site_rate, 
    c(2.427350, 0.0), 
    tolerance = 1e-6
  )
})

test_that("accrue_week: trial structure correct (no capping)", {
  expect_equal(
    aw_out_ls[[2]]@treatment_arm_struct,
    matrix(c(rep(c(TRUE, FALSE, TRUE), each = 2)), ncol = 2)
  )
})


# Testing site_sums()

test_that("site_sums: site totals are correct", {
  expect_equal(
    as.vector(site_sums(aw_out_ls[[1]])),
    c(4, 0)
  )
})


# Testing apply_site_cap()

## Setup

site_cap_accrual_ar <- array(
  data = as.integer(
    c(
      rep(4, 4), 3, rep(0, 7),
      rep(0, 4), 4, rep(0, 7),
      rep(4, 5), rep(0, 7),
      rep(2, 4), 1, rep(0, 7),
      rep(0, 4), 2, rep(0, 7),
      rep(2, 4), 1, rep(0, 7)
    )
  ),
  dim = c(12, 3, 2),
  dimnames = list(
    Weeks = NULL, 
    Arms = c("T1", "T2", "Control"), 
    Centres = paste("Centre", 1:2)
  )
)

fixed_acc_obj@accrual <- site_cap_accrual_ar
fixed_acc_obj@week <- as.integer(5)

fixed_acc_obj <- apply_site_cap(fixed_acc_obj)

test_that("apply_site_cap: Sites with too much accrual are capped.", {
  expect_equal(
    sum(colSums(fixed_acc_obj@accrual[, , 1])),
    fixed_acc_obj@site_cap[1]
  )
})

test_that("apply_site_cap: Capping happens on most recent week only.", {
  expect_equal(
    fixed_acc_obj@accrual[1:4, , 1],
    site_cap_accrual_ar[1:4, , 1]
  )
})

test_that("apply_site_cap: Capping never increases accrual to an arm", {
  expect_true(all(
    fixed_acc_obj@accrual[5, , 1] <= site_cap_accrual_ar[5, , 1]
  ))
})

test_that("apply_site_cap: Sites under the cap are not capped.", {
  expect_equal(
    fixed_acc_obj@accrual[, , 2],
    site_cap_accrual_ar[, , 2]
  )
})

# Testing apply_arm_cap()

fixed_acc_obj@target_arm_size <- as.integer(16)
fixed_acc_obj@target_interim <- as.integer(8)

struct_obj <- trial_structure(
  props_df = data.frame(
    category = c("B1", "B2"),
    region_1 = c(0.54, 0.46)
  ),
  arms_ls = list(T1 = 1:2, T2 = 1:2),
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 2),
    mean_rate = c(8, 4),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 5)
  ),
  precision = NULL,
  shared_control = TRUE,
  control_ratio = c(1, 1),
  fixed_region_prevalences = TRUE
)

#fixed_acc_obj <- apply_arm_cap(fixed_acc_obj)

#print(fixed_acc_obj@accrual)

# Testing get_weeks()

weeks <- get_weeks(c(3, 5))

test_that("get_weeks: produces integer output", {
  checkmate::expect_integer(
    weeks,
    lower = 0,
    any.missing = FALSE,
    len = 2,
    null.ok = FALSE
  )
})

test_that("get_weeks: results are correct", {
  expect_equal(
    weeks,
    c(12, 20)
  )
})

test_that("get_weeks: does not work with negative months", {
  expect_error(get_weeks(-1))
})

test_that("get_weeks: does not work with missing data", {
  expect_error(get_weeks(c(1, NA)))
})

test_that("get_weeks: does not work with null data", {
  expect_error(get_weeks(c(NULL)))
})


