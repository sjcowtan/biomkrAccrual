### Testing do_dirichlet_draws
sites_in_region <- c(1, 2, 3, 2, 1)
region_prevalence <- matrix(c(0.2, 0.6, 0.7, 0.9, 0.4, 0.1), ncol = 3)
precision <- 10

dirichlet_draws_out <- 
  do_dirichlet_draws(region_prevalence, sites_in_region, precision)

test_that("do_dirichlet_draws returns a matrix", {
  checkmate::expect_matrix(
    dirichlet_draws_out,
    any.missing = FALSE,
    nrows = nrow(region_prevalence),
    ncols = length(sites_in_region), 
    null.ok = FALSE,
  )
})

test_that("do_dirichlet_draws contains valid probabilities", {
  checkmate::expect_numeric(
    dirichlet_draws_out,
    lower = 0,
    upper = 1
  )
})

test_that("Columns of do_dirichlet_draws sum to 1", {
  expect_equal(
    colSums(dirichlet_draws_out),
    rep(1, length(sites_in_region))
  )
})

### Testing get_recruit_arm_prevalence

centres_df <- data.frame(
  site = 1:6,
  start_month = 1,
  mean_rate = 10,
  region = c(1, 1, 2, 1, 2, 2),
  site_cap = 70,
  start_week = 1
)

props_df <- data.frame(
  category = LETTERS[1:2],
  proportion_1 = c(0.2, 0.5),
  proportion_2 = c(0.35, 0.9),
  proportion_3 = c(0.01, 0.6)
)

test_that("At least one site", {
  expect_error(get_recruit_arm_prevalence(
    props_df,
    data.frame(
      site = NULL,
      start_month = NULL,
      mean_rate = NULL,
      region = NULL,
      site_cap = NULL,
      start_week = NULL
    ),  
    10,
    fixed_region_prevalences = FALSE
  ))
})

test_that("Sites must be in regions we have prevalences for", {
  expect_error(get_recruit_arm_prevalence(
    c(sample(1:3, 9, replace = TRUE), ncol(props) + 1), 
    props_df, 
    data.frame(
      site = 1:6,
      start_month = 1,
      mean_rate = 10,
      region = c(1, 1, 2, 1, 9, 2),
      site_cap = 70,
      start_week = 1
    ),
    10,
    fixed_region_prevalences = FALSE
  ))
})

test_that("Proportions must not be missing", {
  expect_error(get_recruit_arm_prevalence(
    centres_df,
    data.frame(
      proportion_1 = c(0.2, 0.5),
      proportion_2 = c(0.35, NA),
      proportion_3 = c(0.01, 0.6)
    ),
    10,
    fixed_region_prevalences = FALSE
  ))
})

test_that("Proportions must be valid", {
  expect_error(get_recruit_arm_prevalence(
    data.frame(
      proportion_1 = c(0.2, 0.5),
      proportion_2 = c(0.35, 1.6),
      proportion_3 = c(0.01, 0.6)
    ),
    centres_df,
    10,
    fixed_region_prevalences = FALSE
  ))
})

test_that("Precision must be positive > 0", {
  expect_error(get_recruit_arm_prevalence(props_df, centres_df, 0, TRUE))
})

recruit_arm_prevalence_out <- 
  get_recruit_arm_prevalence(props_df, centres_df, 10, FALSE)

test_that("get_recruit_arm_prevalence returns a matrix", {
  checkmate::expect_matrix(
    recruit_arm_prevalence_out,
    any.missing = FALSE,
    nrows = nrow(props_df),
    ncols = nrow(centres_df),
    null.ok = FALSE
  )
})

test_that("get_recruit_arm_prevalence matrix contains valid probabilities", {
  checkmate::expect_numeric(
    recruit_arm_prevalence_out,
    lower = 0,
    upper = 1
  )
})

test_that("Prevalence sets sum to 1", {
  expect_equal(
    colSums(recruit_arm_prevalence_out), 
    rep(1, nrow(centres_df))
  )
})


### Testing constructor trial_structure()

ts_obj <- trial_structure(
  props_df <- data.frame(
    category = LETTERS[1:2],
    A = c(0.2, 0.8),
    B = c(0.35, 0.65)
  ),
  arms_ls = list(
    T1 = as.integer(1),
    T2 = as.integer(2)
  ),
  centres_df = data.frame(
    site = 1:2,
    start_month = c(1, 5),
    mean_rate = c(10, 18),
    region = c(1, 1),
    site_cap = c(40, 20),
    start_week = c(1, 20)
  ),
  precision = 10,
  shared_control = TRUE,
  fixed_region_prevalences = FALSE
)

test_that(paste(
  "Constructor for trial_structure produces an object of classes",
  "S7_object and biomkrAccrual::trial_structure"
), {
  checkmate::expect_class(ts_obj, c(
    "biomkrAccrual::trial_structure",
    "S7_object"
  ))
})


### Testing is.trial_structure

test_that("is.trial_structure correctly identifies class", {
  expect_true(is.trial_structure(ts_obj))
})

## Testing arm IDs

test_that("treatment_arm_ids has 2 ids, 1 and 2", {
  expect_equal(
    ts_obj@treatment_arm_ids, 
    list(T1 = as.integer(c(1)), T2 = as.integer(c(2)))
  )
})

### Testing arm removal

ts_minusarm_obj <- remove_treat_arms(ts_obj, as.integer(2))

test_that("remove_treat_arms correctly removes arm 2", {
  expect_equal(
    ts_minusarm_obj@treatment_arm_ids, 
    list(T1 = as.integer(c(1)), T2 = c(NA_integer_))
  )
})

test_that("remove_treat_arms changes treatment_arm_struct", {
  expect_false(isTRUE(all.equal(
    ts_minusarm_obj@treatment_arm_struct,
    ts_obj@treatment_arm_struct
  )))
})

test_that("treatment_arm_struct column 2 is now FALSE", {
  expect_false(all(ts_minusarm_obj@treatment_arm_struct[, 2]))
})

test_that("recruit_arm_prevalence has not changed for arm 1", {
  expect_equal(
    ts_minusarm_obj@recruit_arm_prevalence[1, ],
    ts_obj@recruit_arm_prevalence[1, ]
  )
})

test_that("recruit_arm_prevalence has changed for arm 2", {
  expect_equal(
    ts_minusarm_obj@recruit_arm_prevalence[2, ],
    rep(0, ncol(ts_obj@recruit_arm_prevalence)),
    tolerance = 1e-6
  )
})

### Testing get_array_prevalence

arm_struct_mx <- matrix(
  c(F, T, T, F, T, F, T, F, T, F, T, F),
  nrow = 3
)

rec_arm_prev <- matrix(
  c(
    0.1, 0.2, 0.7,
    0.3, 0.1, 0.6
  ),
  nrow = 3
)

control_ratio <- c(0.75, 0.25)

arm_prev_ar <- get_array_prevalence(
  arm_struct_mx, 
  rec_arm_prev, 
  shared_control = FALSE,
  control_ratio
)

test_that("Dimensions of prevalence array are correct", {
  expect_equal(
    dim(arm_prev_ar), 
    c(
      nrow(arm_struct_mx), 
      ncol(arm_struct_mx) * 2, 
      ncol(rec_arm_prev)
    )
  )
})

test_that("Total recruitment percentages are correct for centre 2", {
  expect_equal(
    rowSums(arm_prev_ar[, , 2]),
    rec_arm_prev[, 2], 
    tolerance = 1e-6
  )
})

test_that("Control ratio is correct for centre 2", {
  expect_equal(
    arm_prev_ar[, seq_len(ncol(arm_struct_mx)), 2] / control_ratio[1],
    arm_prev_ar[
      , 
      seq(ncol(arm_struct_mx) + 1, length.out = ncol(arm_struct_mx)),
      2
    ] / control_ratio[2],
    tolerance = 1e-6
  )
})