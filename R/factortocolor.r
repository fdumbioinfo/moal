#' @title change a factor levels with color
#' @description change a factor levels with color
#' @param factor factor
#' @details use an intern define color palette
#' @return character of color value
#' @examples
#' # not run
#' # factortocolor(factor)
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach %do%
#' @noRd
factortocolor <- function( factor = NULL )
{
  factor %>% levels -> levels
  factor %>% as.character -> factorcol
  i = 1
  foreach(i=1:length(levels)) %do% { replace( factorcol , factor == levels[i] , palette0[i] ) -> factorcol }
  factorcol %>% return()
}