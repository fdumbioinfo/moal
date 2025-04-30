#' @title compute confident interval
#' @param y numeric vector
#' @param conf.int numeric
#' @details return se value
#' @return numeric
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@u-psud.fr>
#' @importFrom magrittr %>%
#' @importFrom stats qt
#' @noRd
ci <- function(y , conf.int = 0.05)
{
  qt( 0.95/2 + conf.int, length(y)-1) -> ci0
  se(y) * ci0 -> ci1
  ci1 %>% return()
}