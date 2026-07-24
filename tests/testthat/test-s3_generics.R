# Testing S3 generics

#### Separate control arms

unshared_acc_obj <- accrual(
  treatment_arm_ids = list(T1 = as.integer(1), T2 = as.integer(2)),
  shared_control = FALSE,
  target_df = data.frame(
    arm = c("T1", "T2"),
    interim = as.integer(c(30, 30)),
    final = as.integer(c(60, 60))
  ),
  target_times = c(6, 12),
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

# Test print.accrual

print_out <- capture.output(print(unshared_acc_obj))

test_that("print.accrual: Output has correct format", {
  checkmate::expect_atomic_vector(
    print_out,
    len = 13,
    any.missing = FALSE,
    names = NULL,
    unique = TRUE
  )
})

test_that("print.accrual: Output has correct headers", {
  expect_equal(
    print_out[1],
    "   T1 T2 Control.T1 Control.T2"
  )
})

test_that("print.accrual: Output has correct row content", {
  expect_equal(
    print_out[[13]],
    "12  0  0          0          0"
  )
})

# Test summary.accrual

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
    start_week = c(1, 5)
  ),
  precision = 10,
  shared_control = FALSE,
  control_ratio = c(1, 1),
  fixed_region_prevalences = FALSE
)

# Add 6 weeks' accrual
aw_out_ls <- accrue_week(
  accrual_obj = unshared_acc_obj, 
  struct_obj = ts_obj
)
aw_out_ls <- accrue_week(
  accrual_obj = aw_out_ls[[1]], 
  struct_obj = ts_obj
)
aw_out_ls <- accrue_week(
  accrual_obj = aw_out_ls[[1]], 
  struct_obj = ts_obj
)
aw_out_ls <- accrue_week(
  accrual_obj = aw_out_ls[[1]], 
  struct_obj = ts_obj
)
aw_out_ls <- accrue_week(
  accrual_obj = aw_out_ls[[1]], 
  struct_obj = ts_obj
)
aw_out_ls <- accrue_week(
  accrual_obj = aw_out_ls[[1]], 
  struct_obj = ts_obj
)


summary_out <- capture.output(summary(aw_out_ls[[1]]))

test_that("summary.accrual: Output has correct format", {
  checkmate::expect_atomic_vector(
    summary_out,
    len = 24,
    any.missing = FALSE,
    names = NULL
  )
})

test_that("summary.accrual: Arm names correct", {
  expect_equal(
    summary_out[2],
    "       T1               T2      Control.T1        Control.T2"
  )
})

test_that("summary.accrual: Means correct for recruitment by arm", {
  expect_equal(
    summary_out[6],
    " Mean   :0.1667   Mean   :0   Mean   :0.08333   Mean   :0   "
  )
})

test_that("summary.accrual: Sites correct", {
  expect_equal(
    summary_out[11],
    "    Centre.1       Centre.2"
  )
})

test_that("summary.accrual: Means correct for recruitment by site", {
  expect_equal(
    summary_out[15],
    " Mean   :0.25   Mean   :0  "
  )
})

test_that("summary.accrual: Arm closure weeks correct", {
  expect_equal(
    summary_out[21],
    "        NA         NA         NA         NA " 
  )
})

test_that("summary.accrual: Arm totals correct", {
  expect_equal(
    summary_out[24],
    "[1] 2 0 1 0"
  )
})

# Test plot.accrual

plot_out <- plot(aw_out_ls[[1]], plot = FALSE)

test_that("plot.accrual: Output is a ggplot object", {
  expect_true(is_ggplot(plot_out))
})

# Test print.trial_structure

print_out <- capture.output(print(ts_obj))

test_that("print.trial_structure: Output has correct format", {
  checkmate::expect_atomic_vector(
    print_out,
    len = 4,
    any.missing = FALSE,
    names = NULL,
    unique = TRUE
  )
})

test_that("print.trial_structure: Headers are correct", {
  expect_equal(
    print_out[1],
    "      T1    T2"
  )
})

test_that("print.trial_structure: Structure is correct", {
  expect_equal(
    print_out[-1],
    c(
      "B1  TRUE FALSE",
      "B2  TRUE  TRUE",
      "B3 FALSE  TRUE"
    )
  )
})

# Test plot.trial_structure

plot_out <- plot(ts_obj)

test_that("plot.trial_structure: Output is a ggplot object", {
  expect_true(is_ggplot(plot_out))
})

# Test summary.trial_structure

summary_out <- summary(ts_obj, quietly = TRUE)

test_that("summary.accrual: Output has correct format", {
  expect_data_frame(
    summary_out,
    any.missing = FALSE,
    null.ok = FALSE,
    nrows = 3,
    ncols = 2,
    row.names = "unique",
    col.names = "unique"
  )
})

test_that("summary.accrual: Column names are correct", {
  expect_equal(
    names(summary_out),
    paste("Site", 1:2)
  )
})

test_that("summary.accrual: Row names are correct", {
  expect_equal(
    rownames(summary_out),
    paste0("B", 1:3)
  )
})

test_that("summary.accrual: Prevalences are correct", {
  expect_equal(
    unname(unlist(summary_out)),
    c(0.3783, 0.3456, 0.2761, 0.3923, 0.3445, 0.2631)
  )
})