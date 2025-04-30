#' @title boxplot for one var
#' @param dat data.frame
#' @param xlab character
#' @param ylab character
#' @param log logical if TRUE data are delog in base 2
#' @return plot
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom plyr ddply
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @noRd
boxplot1 <- function( dat, ylab = "y", xlab = "TREATMENT", log=T )
{
  if(log){ dat[,1] %>% as.numeric %>% "^"(2) -> dat[,1] }
  par(xpd=T, mar=c(5,4,5,6))
  p <- ggboxplot(dat, x = ".", y = "y",
                 color = ".", palette =palette0[1:length(dat[,2] %>% table)],
                 add = "jitter")
  dat[,2] %>% levels %>% rev %>% as.character %>% combn( 2 , simplify = F ) -> comp 
  p + stat_compare_means( comparisons = comp ) -> p
  p + stat_compare_means( label.y = (max( dat %>% "["(1) %>% unlist )+0.5) ) -> p
  p + xlab(xlab) -> p
  p + ylab("Intensity") -> p
  p + ggtitle(ylab) -> p
  p + theme( legend.position = c(0.90,0.9),legend.key.size = unit(0.4,'cm'),
             axis.text.x=element_text(face="plain",color="black",size=7,angle=0),
             legend.title = element_text(size=9)) -> p
  p %>% return()
}
