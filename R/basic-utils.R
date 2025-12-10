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