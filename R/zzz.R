#' Enable usage of <S7_object>@name in package code even for 
#' older versions of R.
#' 
#' @name enable-for-old-R
#' @rawNamespace if (getRversion() < "4.3.0") importFrom("S7", "@")
NULL


#' Ensure methods are registered when methods for generics are implemented in 
#' other packages (S3, S4 & S7). This works like export directives in the 
#' NAMESPACE does for S3 and S4.
#' @param ... Parameters passed to S7::methods_register()
#' 
#' @name register-S7-methods
.onLoad <- function(...) {
  S7::methods_register()
}

