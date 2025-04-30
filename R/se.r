#' @title standart error
#' @description
#' Function to make basic line plot from a numeric serie and factor.
#' @param y numeric vector
#' @details return se value
#' @return numeric
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @noRd
se <- function(y)
{
  sd(y) / sqrt( length(y) ) -> sd0
  sd0 / sqrt(length(y)) -> sd1
  sd1 %>% return()
}