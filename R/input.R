#' Convert centres_df to more useable format
#' @param centres_df Dataframe with columns site, start_month, mean_rate,
#' region and site_cap
#' @return Dataframe with start_week, ordered by start_month and site
#' 
do_clean_centres <- function(centres_df) {

  # Get start weeks
  centres_df$start_week <- get_weeks(centres_df$start_month - 1) + 1

  # Order centres_df by start month and site number
  centres_df <- centres_df[with(centres_df, order(start_week, site)), ]

  return(centres_df)
}