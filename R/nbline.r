#' @title count lines in a file
#' @description count lines in a file
#' @param filepath character
#' @param compression character read file compression type
#' @examples
#' # not run
#' # nbline( filepath )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @noRd
nbline <- function( filepath , compression = "none"  )
{
  ### none
  if( compression == "none" )
  {
    filepath %>% file("r") -> f
    i=0
    while(T)
    {
      readLines(f, n=1 ) -> line
      if( length(line) == 0 ){ break }
      i=i+1
    }
    close(f)
  }
  ### gz
  if( compression == "gz" )
  {
    filepath %>% gzfile("r") -> f
    i=0
    while(T)
    {
      readLines(f, n=1 ) -> line
      if( length(line) == 0 ){ break }
      i=i+1
    }
    close(f)
  }
  ### zip
  if( compression == "zip" )
  {
    filepath %>% unz("r") -> f
    i=0
    while(T)
    {
      readLines(f, n=1 ) -> line
      if( length(line) == 0 ){ break }
      i=i+1
    }
    close(f)
  }
  i=i-1
  return(i)
}