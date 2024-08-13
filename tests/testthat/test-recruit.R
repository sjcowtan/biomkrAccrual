### Testing accrual() constructor

acc_obj <- accrual(
  treatment_arm_ids = list(T1 = as.integer(1), T2 = as.integer(2)),
  shared_control = TRUE,
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 5),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 20)
  ),
  accrual_period = as.integer(36)
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
