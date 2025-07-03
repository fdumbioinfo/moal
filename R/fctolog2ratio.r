#' @title Funtion to transform fold-change into log2 ratio
#' @description Funtion to transform fold-change into log2 ratio
#' @param fc numeric vector with fold-change value
#' @details usefull to format data for volcano plot
#' @return numeric vector with logratio value
#' @examples
#' library(magrittr)
#' c(1.5 , -1.5 ,5, 2 , -10 ) %>% fctolog2ratio
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @noRd
fctolog2ratio <- function( fc )
{
  fc %>% as.numeric %>% "<"(.,0) %>% which -> sel
  if(length(sel)>0){ fc[sel] %>% invratio %>% "*"(.,-1) -> fc[sel] }
  fc %>% log2 %>% return()
}