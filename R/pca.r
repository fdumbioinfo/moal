#' @title Principal Component Analysis
#' @description Principal Component Analysis
#' @param dat matrix numeric
#' @param factor factor
#' @param samplename character
#' @param pc1 numeric
#' @param pc2 numeric
#' @param center logical TRUE
#' @param scale logical TRUE
#' @param title character
#' @param legendtitle character
#' @return pca plot
#' @examples
#' # not run
#' # pca(dat,factor)
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom scales percent
#' @noRd
pca <- function(
  dat = NULL, factor = NULL, samplename = NULL, pc1 = 1, pc2 = 2,
  center = TRUE ,scale = TRUE, title = "PCA", legendtitle = "TREATMENT" )
{
  # data
  if(!is.null(samplename)){ dat %>% setNames(samplename) -> dat }
  dat %>% as.matrix %>% t %>% stats::prcomp(center=center,scale. = scale) -> pca0
  pca0$x %>% data.frame -> pca1
  prop <- pca0$sdev^2/sum(pca0$sdev^2)
  # plot
  if(is.null(factor))
  {
    pca1 %>% dplyr::select(c(pc1,pc2)) %>% setNames(c("PC1","PC2")) -> pca2
    pca2 %>% ggplot(aes(x=.data$PC1,y=.data$PC2)) -> p
    # pca1 %>% ggplot(aes(x=colnames(pca1)[pc1],y=colnames(pca1)[pc2])) -> p
    p + geom_point(size=5) -> p
    p + geom_text(label=colnames(dat),size=2,hjust=0,vjust=3,color="black") -> p
    p + labs(
      x = paste("PC",pc1,"(",scales::percent(prop[pc1]),")",sep=""),
      y = paste("PC",pc2,"(",scales::percent(prop[pc2]),")",sep="")) -> p
    p + guides(color=guide_legend(legendtitle),shape=guide_legend(legendtitle)) -> p
    p + ggtitle(title) -> p
    p + theme(plot.title=element_text(size=12,face="bold",hjust=0.5),legend.position='right') -> p
    p %>% return()
  }else
    {
      # pca1 %>% ggplot(aes_string(x=paste("PC",pc1,sep=""),y=paste("PC",pc2,sep=""))) -> p
      pca1 %>% dplyr::select(c(pc1,pc2)) %>% setNames(c("PC1","PC2")) -> pca2
      pca2 %>% ggplot(aes(x=.data$PC1,y=.data$PC2)) -> p
      p + geom_point(aes(color=factor),size=5) -> p
      p + geom_text(label=colnames(dat),size=2,hjust=0,vjust=3,color="black") -> p
      p + labs(
        x=paste("PC",pc1,"(",scales::percent(prop[pc1]),")",sep=""),
        y=paste("PC",pc2,"(",scales::percent(prop[pc2]),")",sep="")) -> p
      p + guides(color=guide_legend(legendtitle), shape=guide_legend(legendtitle)) -> p #to change legend name
      p + ggtitle(title) -> p
      p + theme(
        plot.title=element_text(size=12,face="bold",hjust=0.5),
        legend.position='right',
        panel.background=element_rect(colour="black",fill= "white"),
        panel.grid.major.x=element_line(colour="black",linewidth=0.1,linetype="dashed"),
        panel.grid.major.y=element_line(colour="black",linewidth=0.1,linetype="dashed")) -> p
      if(is.numeric(factor)){ p + scale_color_gradient(low="azure",high="blue4") -> p }else{
        p + scale_color_manual(values=palette0) -> p }
      p %>% return()
    }
}
