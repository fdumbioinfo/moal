#' @title Boxplot
#' @param dat matrix numeric
#' @param factor factor
#' @param outline logical display outliers FALSE by default
#' @param title character
#' @param legendtitle character
#' @param outlier boolean
#' @param coefiqr numeric
#' @param ggplot logical use graphics library or ggplot FALSE by default
#' @details To make boxplot from matrix.
#' @return plot
#' @examples
#' # not run
#' # mat1 %>% boxplot( factor = sif1$F3 )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom graphics boxplot legend
#' @importFrom utils stack
#' @noRd
boxplot <- function(
  dat, factor, outline = FALSE, title = "Boxplot", legendtitle = "TREATMENT",
  outlier = T, coefiqr = 1.5, ggplot = FALSE )
{
  # data
  factor %>% order -> sel
  dat[ , sel] -> dat
  factor[sel] -> factor
  #
  if(ggplot)
  {
    dat %>% stack %>% 
      data.frame(TREATMENT=rep(as.character(factor),rep(dim(dat)[1],dim(dat)[2]))) -> dat0
    dat0$TREATMENT %>% ordered( levels( factor ) ) -> dat0$TREATMENT
    # plot
    dat0 %>% ggplot(aes(y=.data$values,x=.data$ind,color=.data$TREATMENT)) -> p
    p + ggtitle(title) -> p
    if(!outline)
    {
      p + ylim(min(dat0$values),min(dat0$values) + 
                 coefiqr*(quantile(dat0$values)[4] - quantile(dat0$values)[2])) -> p
      p + geom_boxplot(coef=coefiqr,outlier.shape=NA) -> p
    }else
      { p + geom_boxplot(coef=coefiqr) -> p }
    p + theme(axis.text.x=element_text(face="plain",color="black",size=1,angle=90)) -> p
    p + theme(axis.text.y=element_text(face="plain",color="black",size=5,angle=0)) -> p
    if(!is.null(factor))
    { p + scale_color_manual(values=palette0) -> p }
    p + guides(color=guide_legend(legendtitle)) -> p
  }
dat %>% graphics::boxplot(main=title,outline=outline,las=2,
                          col=factor %>% factortocolor,cex.axis=0.65)
legend("topright",legend=levels(factor),title=legendtitle,
       col=palette0[1:length(levels(factor))],lty=1,lwd=5,cex=0.7)
}