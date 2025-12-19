#' @title Histogram
#' @description function to make histogram plot
#' @param dat matrix numeric
#' @param factor factor
#' @param title character
#' @param bins numeric
#' @param color character
#' @details To make Histogram plot
#' @return no value
#' @examples
#' # not run
#' # histogram(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom utils stack
#' @importFrom rlang .data
#' @noRd
histogram <- function( dat = NULL , factor = "" , title = "histogram" , bins = 30 , color = "black" )
{
  dat %>% utils::stack(.) -> dat0
  dat0 %>% ggplot(aes(x = .data$values)) -> p
  p + ggtitle(title) -> p
  p + geom_histogram(color=color,fill="white",bins=bins) -> p
  p + theme(axis.text.x=element_text(face="plain",color="black",size=7,angle=90)) -> p
  p + theme(axis.text.y=element_text(face="plain",color="black",size=8,angle=0)) -> p
  p %>% return()
 }