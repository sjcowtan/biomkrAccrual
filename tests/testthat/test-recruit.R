### Testing accrual() constructor

acc_obj <- accrual(
  treatment_arm_ids = list(T1 = as.integer(1), T2 = as.integer(2)),
  shared_control = TRUE,
  accrual_period = as.integer(12),
  interim_period = as.integer(6),
  control_ratio = c(1, 1),
  var_lambda = 0.25,
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 5),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 20)
  )
)

test_that(paste(
  "Constructor for accrual produces an object of classes",
  "S7_object and biomkrAccrual::accrual"
), {
  checkmate::expect_class(acc_obj, c(
    "biomkrAccrual::accrual",
    "S7_object"
  ))
})


### Testing is.accrual()

test_that(paste(
  "Constructor for accrual produces an object of classes",
  "S7_object and biomkrAccrual::accrual"
), {
  checkmate::expect_class(acc_obj, c(
    "biomkrAccrual::accrual",
    "S7_object"
  ))
})


### Testing treat_sums()
arr <- array(1:24, 2:4)

test_that("Can sum by treatment a 3-D accrual array", {
  checkmate::expect_integer(
    treat_sums(arr),
    min.len = 1,
    max.len = ncol(arr),
    lower = 0,
    any.missing = FALSE,
    null.ok = FALSE
  )
})

test_that("Output correct from treat_sums for valid array", {
  expect_identical(
    treat_sums(arr),
    as.integer(c(84, 100, 116))
  )
})

## Testing treat_sums on an accrual object

test_that("Output correct from treat_sums for valid accrual object", {
  expect_identical(
    treat_sums(acc_obj),
    as.integer(c(0, 0, 0))
  )
})

## Testing capping
uncapped <- c(1, 1, 2, 5, 8)
capped <- do_choose_cap(uncapped, 3)

test_that("do_choose_cap caps population correctly", {
  expect_identical(
    length(capped),
    as.integer(3)
  )
})

# Need to be able to compare frequencies for elements in capped only
table_uncapped <- as.data.frame(table(uncapped))
table_capped <- as.data.frame(table(capped))

test_that("do_choose_caps does not produce extra repeats", {
  expect_true(
    all(
      table_uncapped[
        which(table_uncapped$uncapped %in% table_capped$capped),
      ]$Freq >= table_capped$Freq
    )
  )
})
