#' @title summarise data with one factor
#' @param dat data.frame
#' @param model character
#' @param log logical
#' @details return ci mean length sd of data.frame variable by factor
#' @return data.frame
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom stats aggregate as.formula sd
#' @noRd
dsum <- function( dat, model, log=T )
{
  if( log ){ dat[,1] %>% "^"(2,.) -> dat[,1] }
  dat[,-1] %>% table %>% as.data.frame -> datN
  paste( names(dat)[1] , "~" , model ) %>% as.formula(.) -> model1
  # mean
  stats::aggregate( y ~ ., data=dat, FUN=mean ) -> datMean
  # sd
  stats::aggregate( y ~ ., data=dat, FUN=sd ) -> datSd
  # se
  stats::aggregate( y ~ ., data=dat, FUN=se ) -> datSe
  # ci
  stats::aggregate( y ~ ., data=dat, FUN=ci ) -> datCi
  #
  data.frame( datMean, datSd[,length(datSd)], datSe[,length(datSe)],
              datCi[,length(datCi)], datN[,length(datN)]) -> dat1
  #
  c( dat %>% names %>% "["(-1),"mean","sd","se","ci","n" ) -> Names2
  names(dat1) <- Names2
  dat1 %>% return()
}