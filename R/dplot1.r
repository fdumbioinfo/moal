#' @title line plot
#' @param dat data.frame
#' @param title character
#' @param legendtitle character
#' @param conf.interval numeric
#' @param log logical if TRUE data are delog in base 2
#' @param ylab character
#' @details return ci mean length sd of data.frame variable by factor
#' @return plot
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@u-psud.fr>
#' @importFrom plyr ddply
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @noRd
dplot1 <- function(
  dat, title="plot", legendtitle="legend", conf.interval=.95, log=T, ylab="y" )
{
  names( dat ) -> Names
  dat %>% dsum( model = Names[2] , log ) -> dat1
  names( dat1 )[1] <- "f"
  #
  dat1 %>% ggplot( aes( x=.data$f, y=.data$mean, color=.data$f, group=1) ) -> p
  p + geom_errorbar( aes( ymin=mean-se, ymax=mean+se ), colour="black", width=.1, position=position_dodge(0.1) ) -> p
  p + geom_line( position = position_dodge( 1 ) ) -> p
  p + geom_point( position=position_dodge(0.1), size=3, shape=21, fill="white" ) -> p
  p + ylab( ylab ) -> p
  p + xlab( title ) -> p
  p + scale_color_manual( values=palette0, name=Names[2] ) -> p
  p + ggtitle( title ) -> p
  p + expand_limits( y=0 ) -> p
  p + theme_bw( ) -> p
  p + theme( legend.justification=c(1,0), legend.position= c(1,0) ) -> p
  p %>% return()
}