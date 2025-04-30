#' @title Line plot
#' @description make line plot from a numeric serie and factor
#' @param dat data.frame
#' @param title character
#' @param legendtitle character
#' @param conf.interval numeric
#' @param log logical if TRUE data are delog in base 2
#' @param ylab character
#' @details return ci mean length sd of data.frame variable by factor
#' @return plot
#' @author Florent Dumont <florent.dumont@universite-paris-sacly.fr>
#' @importFrom plyr ddply
#' @importFrom magrittr %>%
#' @noRd
dplot <- function( dat, title="plot", legendtitle="legend", conf.interval=.95 , log=T, ylab="y" )
{
  if( length(dat) > 2 )
  {
    dplot2( dat=dat, title=title, legendtitle=legendtitle,
            conf.interval = conf.interval, log = log, ylab = ylab )
  }else{
    dplot1( dat=dat, title=title, legendtitle=legendtitle,
            conf.interval=conf.interval, log=log, ylab=ylab ) }
}
