#' @title line plot
#' @description
#' Function to make basic line plot from a numeric serie with One factor
#' @param dat data.frame
#' @param title character
#' @param legendtitle character
#' @param conf.interval numeric
#' @param log logical if TRUE data are delog in base 2
#' @param ylab character
#' @details return ci mean length sd of data.frame variable by factor
#' .model exemple = "GENOTYPE+TREATMENT"
#' @return plot
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@u-psud.fr>
#' @importFrom plyr ddply
#' @importFrom magrittr %>%
#' @noRd
dplot2 <- function(
  dat,
  title = "plot" ,
  legendtitle = "legend",
  conf.interval = .95 ,
  log = T,
  ylab = "y" )
{
  names(dat) -> Names
  paste( Names[2],"*",Names[3] ) -> Model
  dat %>% dsum( model = Model, log = log ) -> dat1
  names(dat1)[c(1,2)] <- c("f1","f2")
  dat1 %>%
    ggplot(
      aes(
        x = .data$f1,
        y = .data$mean,
        color = .data$f2,
        group = .data$f2) ) -> p

  p + geom_errorbar(
    aes( ymin = mean-se , ymax = mean+se ),
    colour = "black", width = .1,
    position = position_dodge( 0.1 ) ) -> p

  p + geom_line( position = position_dodge( 0.1 ) ) -> p
  p + geom_point(
    position = position_dodge(0.1),
    size = 3 , shape = 21 , fill = "white" ) -> p
  p + ylab( ylab ) -> p
  p + xlab( Names[2] ) -> p
  p + scale_color_manual(
    values = palette0,
    name = Names[3] ) -> p
  p + ggtitle( title ) -> p
  p + expand_limits( y = 0 ) -> p
  p + theme_bw() -> p
  p + theme(
    legend.justification = c(1,0),
    legend.position = c(1,0) ) -> p
  p %>% return()
}