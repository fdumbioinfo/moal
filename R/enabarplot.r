#' @title ena barplot
#' @param dat data.frame enrichment analysis results with 3 column
#' @param datnes data.frame enrichment analysis results with 3 column
#' @param labsize numeric feature size
#' @param addratioena logical if TRUE add overlap and geneset size on enrichment barplot
#' @param addenarankbarplot logical if TRUE add ena barplot ranked by NES score
#' @param title character
#' @param subtitle character
#' @return plot
#' @examples
#' # not run
#' # enabarplot( dat )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats reorder
#' @noRd
enabarplot <- function( dat = NULL,datnes = NULL, labsize = 7, addratioena = TRUE, addenarankbarplot = TRUE, title = "ena", subtitle = "")
{
  # barplot
  
  
  
  
  
  # fgseapval1plot3 %>% colnames
  dat %>% ggplot( aes(x=Log10Pval,y=.data$Name,fill=.data$NES)) -> p
  p + geom_bar(stat="identity",width = 0.95) -> p
  p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
  # p + theme_minimal() -> p
  # p + theme_light() -> p
  # p + theme_linedraw() -> p
  p + theme_bw() -> p
  p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 5),
            axis.text.y = element_text(angle = 0, vjust = 0.4, hjust=1,size = 6)) -> p
  p + labs(x="log10(p-value)",y="Genesets") -> p
  title %>% paste(" (Log10p ranking)",sep="") -> Title0
  # p + ggtitle(Title0) -> p
  p + ggtitle(Title0,subtitle=subtitle) -> p
  #
  if(addratioena){p + geom_text(data=dat,aes(label=.data$Overlap),hjust=1.1,size=2) -> p}
  p + theme(plot.title=element_text(size=10,hjust=0.5),
            plot.subtitle=element_text(size=8,hjust=0.5),
            axis.title.x=element_text(size=6,face="bold"),
            axis.text.x = element_text(size=6,face="bold"),
            axis.title.y=element_text(size=6,face="bold"),
            axis.text.y=element_text(size=8,face="bold"),
            legend.position="bottom",
            legend.title=element_text(size = 5),
            legend.text=element_text(size = 5)) -> p
  p + geom_vline(xintercept=-log10(0.05),color="darkgoldenrod1",size=0.4,linetype="dashed") -> p
  p -> p1
  if(addenarankbarplot)
  {
    if(!is.null(datnes)){ datnes -> dat }
    dat %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,abs(.data$NES))) -> dat
    # dat$OverlapSize %>% paste("/",dat$GenesetSize) %>% 
    #   lapply(paste0,collapse="") %>% unlist -> Overlap
    # fgseapval1plot3 %>% data.frame(Overlap) -> fgseapval1plot3
    dat %>% ggplot( aes(x=.data$Log10Pval,y=.data$Name,fill=.data$NES)) -> p
    p + geom_bar(stat="identity") -> p
    p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
    p + theme_bw() -> p
    if(addratioena){p + geom_text(data = dat, aes(label = .data$Overlap),hjust=1.1,size = 2) -> p}
    p + theme(plot.title=element_text(size=10,hjust=0.5),
              plot.subtitle=element_text(size=8,hjust=0.5),
              axis.title.x=element_text(size=6,face="bold"),
              axis.text.x=element_text(size=6,face="bold"),
              axis.title.y=element_text(size=6,,face="bold"),
              axis.text.y=element_text(size=8,face="bold"),
              legend.position="bottom",
              legend.title=element_text(size = 5), 
              legend.text=element_text(size = 5)) -> p
    p + labs(x="log10(p-value)",y="Genesets") -> p
    title %>% paste(" (NES ranking)",sep="") -> Title0
    # p + ggtitle(Title0) -> p
    p + ggtitle(Title0,subtitle=subtitle) -> p
    p + geom_vline(xintercept=-log10(0.05),color="darkgoldenrod1",size=0.4,linetype="dashed") -> p
    p -> p2
    # gridExtra::grid.arrange( p1, p2, nrow = 1 ) -> p
    # gridExtra::grid.arrange(p1,p2,heights = c(5, 5),nrow=2,) -> p
    # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),nrow=2) -> p
    # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),nrow=2) -> p
    # gridExtra::grid.arrange(p1,p2,heights = c(10, 10),ncol=2) -> p
    # gridExtra::grid.arrange(p1,p2,heights = c(5,5),ncol=2) -> p
    gridExtra::grid.arrange(p1,p2,ncol=2) -> p
    # gridExtra::grid.arrange(p1,p2,heights=c(10, 10)) -> p
    # ggpubr::ggarrange(p1, p2, widths = c(5,5)) -> p
    # plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
  }
  p %>% return(.)
}
