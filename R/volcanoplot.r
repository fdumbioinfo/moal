#' @title Volcanoplot
#' @description Do volcanoplot
#' @param dat data.frame table with 4 columns (see details)
#' @param pval numeric p-value threshold
#' @param fc numeric fold-change threshold
#' @param topgenename logical display gene label TRUE by default
#' @param topgenenamen numeric increase number of gene label
#' @param genenamelist character vector of gene list to label
#' @param genenamesize numeric label size for gene name
#' @param title character
#' @details
#' dat parameter must have 4 columns: rowID , p_AvsB , fc_AvsF and Symbol.
#' @return no returned value
#' @examples
#' # not run
#' # data.frame(rowID,p_AvsB,fc_AvsB,Symbol) -> dat
#' # volcanoplot(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr slice
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @export
volcanoplot <- function(
    dat = NULL, pval = 0.05, fc = 1.5, topgenename = TRUE, topgenenamen = 5,
    genenamelist = NULL, genenamesize = 2, title = "Volcanoplot" )
{
  dat -> Dat0
  dat[,2] %>% log10 %>% "*"(.,-1) -> Dat0[,2]
  dat[,3] %>% fctolog2ratio -> Dat0[,3]
  #
  ggplot() -> p
  # black color for all point
  p + geom_point(data=Dat0,aes(x=.data[[colnames(Dat0)[3]]],y=.data[[colnames(Dat0)[2]]]),color="black") -> p
  # red color for up point
  Dat0 %>% dplyr::filter(.[[colnames(Dat0)[2]]] %>% ">"(-log10(pval)) & .[[colnames(Dat0)[3]]] > log2(fc)) -> u
  if(nrow(u)>0){ p + geom_point(data=u,aes(x=.data[[colnames(u)[3]]],y=.data[[colnames(u)[2]]]),color="red") -> p }
  # green color for down point
  Dat0 %>% dplyr::filter(.[[colnames(Dat0)[2]]] %>% ">"(-log10(pval)) & .[[colnames(Dat0)[3]]] < -log2(fc)) -> d
  if(nrow(d)>0){ p + geom_point(data=d,aes(x=.data[[colnames(d)[3]]],y=.data[[colnames(d)[2]]]),color="green") -> p }
  # -----
  # top genename label
  # -----
  if(topgenename & is.null(genenamelist))
  {
    # label for up point
    if(nrow(u) > 0)
    { 
      u %>% dplyr::arrange(-.[[colnames(u)[2]]]) %>% dplyr::slice(1:topgenenamen) -> up
      u %>% dplyr::arrange(-.[[colnames(u)[3]]]) %>% dplyr::slice(1:topgenenamen) -> ufc
      rbind(up,ufc) %>% unique -> uname
      p + ggrepel::geom_label_repel(data=uname,aes(x=.data[[colnames(uname)[3]]],y=.data[[colnames(uname)[2]]],
                                          label=.data[[colnames(uname)[4]]]),
                           size=genenamesize,color="black") -> p
    }
    # label for down point
    if(nrow(d) > 0)
    { 
      d %>% dplyr::arrange(-.[[colnames(d)[2]]]) %>% dplyr::slice(1:topgenenamen) -> dp
      d %>% dplyr::arrange(-.[[colnames(d)[3]]]) %>% dplyr::slice(1:topgenenamen) -> dfc
      rbind(dp,dfc) %>% unique -> dname
      p + ggrepel::geom_label_repel(data=dname,aes(x=.data[[colnames(dname)[3]]],y=.data[[colnames(dname)[2]]],
                                          label=.data[[colnames(dname)[4]]]),
                           size=genenamesize,color="black") -> p
    }
  }
  # -----
  # genename list label
  # -----
  if(!is.null(genenamelist))
  {
    genenamelist %>% as.character %>% unique %>% data.frame(Symbol=.) %>% 
      dplyr::inner_join(Dat0,by="Symbol") %>% dplyr::select(c(2,3,4,1)) -> genename
    #
    if(nrow(genename) > 0)
    {
      p + geom_point(data=genename,aes(x=.data[[colnames(genename)[3]]],y=.data[[colnames(d)[2]]]),color="orange") -> p 
      # red color for up point
      Dat0 %>% dplyr::filter(.[[colnames(Dat0)[2]]] %>% ">"(-log10(pval)) & .[[colnames(Dat0)[3]]] > log2(fc)) -> u
      if(nrow(u)>0){ p + geom_point(data=u,aes(x=.data[[colnames(u)[3]]],y=.data[[colnames(u)[2]]]),color="red") -> p }
      # green color for down point
      Dat0 %>% dplyr::filter(.[[colnames(Dat0)[2]]] %>% ">"(-log10(pval)) & .[[colnames(Dat0)[3]]] < -log2(fc)) -> d
      if(nrow(d)>0){ p + geom_point(data=d,aes(x=.data[[colnames(d)[3]]],y=.data[[colnames(d)[2]]]),color="green") -> p }
      #
      p + ggrepel::geom_label_repel(data=genename,
                                    aes(x=.data[[colnames(genename)[3]]],y=.data[[colnames(genename)[2]]],label=.data[[colnames(genename)[4]]]),
                                    size=genenamesize,color="black",max.overlaps=nrow(genename)) -> p
      
    }
  }
  #
  p + theme_bw() -> p
  p + theme(plot.title=element_text(size=10,face="bold",hjust=0.5),legend.position="none") -> p
  p + geom_hline(yintercept=-log10(pval),linetype="solid",colour="orange",linewidth=0.3) -> p
  p + geom_vline(xintercept=log2(fc),linetype="solid",colour="orange",linewidth=0.3) -> p
  p + geom_vline(xintercept=-log2(fc),linetype="solid",colour="orange",linewidth=0.3) -> p
  p + labs(x=paste("log2ratio(",fc,")",sep=""),y=paste("-log10(",pval,")",sep="")) -> p
  p + ggtitle(paste(title, sep="")) -> p
  p %>% return()
}
