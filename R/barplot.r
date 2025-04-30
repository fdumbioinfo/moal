#' @title Barplot
#' @param dat data.frame with 2 column
#' @param title character
#' @param titlesize numeric
#' @param xlab character
#' @param ylab character
#' @param xsize numeric
#' @param ysize numeric
#' @details make barplot
#' @return plot
#' @examples
#' # not run
#' # barplotGO(goEnrichment)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom graphics boxplot legend
#' @noRd
barplot <- function(
    dat,
    title = "GO Barplot",
    titlesize = 10,
    xsize = 5,
    ysize = 5,
    xlab = "Biological process",
    ylab = "Enrichment score")
{
  dat %>% setNames(c("Term","values")) -> dat
  ggplot( dat , aes( x = dat$Term, y = dat$values )) +
    stat_summary(geom = "bar", fun = mean, position = "dodge") +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    scale_y_continuous(breaks = round(seq(0, max(dat$values), by = 2), 1)) +
    theme_bw(base_size=20) +
    theme(
      legend.position='none',
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=titlesize, face="bold", vjust=1),
      axis.text.x=element_text(angle=0, size=xsize, face="bold", hjust=1.10),
      axis.text.y=element_text(angle=0, size=ysize, face="bold", vjust=0.5),
      axis.title=element_text(size=5, face="bold"),
      legend.key=element_blank(),     #removes the border
      legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
      legend.text=element_text(size=5),  #Text size
      title=element_text(size=5) ) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip() -> p
  p %>% return()
}
