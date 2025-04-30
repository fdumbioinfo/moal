.onLoad <- function (libname, pkgname)
{
  # prevent note about global variables in devtools::check()
  utils::globalVariables (".")
  invisible ()
}