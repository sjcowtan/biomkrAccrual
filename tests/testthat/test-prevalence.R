### Testing do_dirichlet_draws
sites_in_region <- c(1, 2, 3, 2, 1)
region_prevalence <- matrix(c(0.2, 0.6, 0.7, 0.9, 0.4, 0.1), ncol = 3)
precision <- 10

dirichlet_draws_out <- do_dirichlet_draws(region_prevalence, sites_in_region, precision)

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
