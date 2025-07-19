#' @title Normalization function
#' @description normalization and log2
#' @param dat data.frame
#' @param method character apply quantile normalization by default see details
#' @param log logical apply log base 2
#' @return data.frame
#' @details
#' for see limma normalizeBetweenArrays method
#' @examples
#' # not run
#' # norm(dt)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom limma normalizeBetweenArrays
#' @export
norm <- function( dat , method = NULL , log = TRUE )
{
  if(is.null(method)){method <- "quantile"}
  dat %>% normalizeBetweenArrays( method=method ) %>% as.data.frame -> dat
  if( log ){ dat %>% "+"(.,1.01) %>% log2 -> dat }
  dat %>% return()
}