#' Create directory if it doesn't already exist, or
#' issue error if it does exist but the minimum access
#' requirements are not met, to avoid using directories
#' created for another purpose.
#' 
#' @param file_path Location of directory to check/create
#' 
#' @importFrom checkmate test_directory_exists assertTRUE
#' @importFrom R.utils fileAccess
#' 

makeifnot_dir <- function(file_path) {
  # Set up directory if does not already exist
  if (dir.exists(file.path(
    file_path
  ))) {
    # Can we write to it?
    ## checkmate::assert_access() relies on base::file.access()
    ## which is unreliable on ubuntu and network drives
    ## c.f. https://github.com/mllg/checkmate/issues/267
    checkmate::assertTRUE(
      all(
        sapply(
          # Exists, read, write, execute
          c(0, 1, 2, 4), 
          function(p)  R.utils::fileAccess(file_path, mode = p)
        ) == 0
      )
    ) 
  } else {
    # Wrap me in a tryCatch
    tryCatch(
      dir.create(file_path, mode = "0777"),
      # Convert warning to error
      warning = function(w) rlang::abort(w$message)
    )
  }
}


#' Expand target_df to include control arm(s)
#' 
#' @param target_df Dataframe of recruitment targets, as read from target file.
#' @param shared_control One shared control arm or one per experimental arm 
#' (logical, defaults to TRUE)?
#' #' @param control_ratio Ratio of patient allocation to treatment arm
#' versus control for all active arms; defaults to c(1, 1).
#' 
expand_targets <- function(
  target_df, 
  shared_control = TRUE, 
  control_ratio = c(1, 1)
) {
  # Want control ratio in form 1:x
  control_ratio <- control_ratio / control_ratio[1]

  if (shared_control) {
    control_df <- as.data.frame(append(
      list(arm = "Control"),
      round(colSums(target_df[-1]) * control_ratio[2], 0)
    ))
  } else {
    control_df <- as.data.frame(append(
      list(arm = paste(target_df$arm, "Control")),
      round(target_df[-1] * control_ratio[2], 0)
    ))
  }

  cat(names(target_df))
  cat(names(control_df))

  rbind(target_df, control_df)
}


#' Convert targets dataframe to long format
#' 
#' @param target_df Dataframe of targets in target file or extended format
#' 
targets_tolong <- function(target_df) {
  target_long_df <- reshape(
    target_df, 
    direction = "long", 
    varying = list(names(target_df)[2:ncol(target_df)]), 
    idvar = "arm", 
    v.names = "value", 
    timevar = "target", 
    times = names(target_df)[-1]
  )

  target_long_df$target <- as.factor(target_long_df$target)

  return(target_long_df)
}