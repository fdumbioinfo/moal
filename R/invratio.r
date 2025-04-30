#' @title invratio
#' @description Function that reverse ratio values
#' @param ratio numeric vector
#' @details
#' a/b become b/a
#' @return numeric vector of inversed ratio
#' @examples
#' # not run
#' # invratioc(0.05 , 0.5 , 2 )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @noRd
invratio <- function( ratio  )
{ 1/ratio %>% return() }