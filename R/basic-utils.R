#' Create directory if it doesn't already exist, or
#' issue error if it does exist but the minimum access
#' requirements are not met, to avoid using directories
#' created for another purpose.
#' 
#' @param file_path Location of directory to check/create
#' 
#' @importFrom checkmate test_directory_exists assert_access
#' 

makeifnot_dir <- function(file_path) {
  # Set up directory if does not already exist
  if (dir.exists(file.path(
    file_path
  ))) {
    # Can we write to it?
    checkmate::assert_access(file.path(
      file_path, access = "rwx"
    )) 
  } else {
    # Wrap me in a tryCatch
    tryCatch(
      dir.create(file_path, mode = "0777"),
      # Convert warning to error
      warning = function(w) rlang::abort(w$message)
    )
  }
}