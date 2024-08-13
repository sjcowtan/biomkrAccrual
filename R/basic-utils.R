#' Create directory if it doesn't already exist, or
#' issue error if it does exist but the minimum access
#' requirements are not met, to avoid using directories
#' created for another purpose.
#' 
#' @param file_path Location of directory to check/create
#' @param min_access Minimum access permissions ("rwx", "rx" etc.)
#' 
#' @importFrom checkmate test_directory_exists assert_access
#' 

makeifnot_dir <- function(file_path, min_access) {
  # Set up directory if does not already exist
  if (checkmate::test_directory_exists(file.path(
    file_path
  ))) {
    # Can we write to it?
    checkmate::assert_access(file.path(
      file_path, access = min_access
    )) 
  } else {
    # Wrap me in a tryCatch
    tryCatch(
      dir.create(file_path),
      # Convert warning to error
      warning = function(w) rlang::abort(w$message)
    )
  }
}