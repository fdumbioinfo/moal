#' @title Volcanoplot
#' @description Do volcanoplot
#' @param dat data.frame table with 4 columns (see details)
#' @param pval numeric p-value threshold
#' @param fc numeric fold-change threshold
#' @param dogenename logical displya GeneName or not. TRUE by default
#' @param GeneNameN logical the number of gene to be displayed. 50 by default
#' @param GeneNameList character vector of gene Symbol. see details.
#' @param GeneNameSize numeric
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
#' @export
volcanoplot <- function(
    dat = NULL, pval = 0.05, fc = 1.5, dogenename = TRUE,
    GeneNameN = 5 , GeneNameList = NULL, GeneNameSize = 2, title = "Volcanoplot" )
{
  dat[,2] %>% log10 %>% "*"(.,-1) -> dat[,2]
  dat[,3] %>% fctolog2ratio -> dat[,3]
  # up
  dat %>% dplyr::slice(which(dat[,2] %>% ">"(.,-log10(pval)) & dat[,3] > log2(fc))) -> Up
  # down
  dat %>% dplyr::slice(which(dat[,2] %>% ">"(.,-log10(pval)) & dat[,3] < -log2(fc))) -> Down
  #
  ggplot() -> p
  p + geom_point(aes(x=dat[,3], y=dat[,2]), colour="black") -> p
  # up
  p + geom_point(aes(x=Up[,3], y=Up[,2]), colour="red") -> p
  # down
  p + geom_point(aes(x=Down[,3], y=Down[,2]), colour="green4") -> p
  if(dogenename)
  {
    if(is.null(GeneNameList))
    {
      # display GeneName up
      # up
      if(nrow(Up)>0)
      { 
        Up %>% dplyr::arrange(-.[[colnames(Up)[2]]]) %>% dplyr::slice(1:GeneNameN) -> GeneNameUpp
        Up %>% dplyr::arrange(-.[[colnames(Up)[3]]]) %>% dplyr::slice(1:GeneNameN) -> GeneNameUpfc
        rbind(GeneNameUpp,GeneNameUpfc) %>% unique -> GeneNameUp
        p + geom_text(aes(x=GeneNameUp[,3],y=GeneNameUp[,2]),
                      label=GeneNameUp[,4],size=GeneNameSize,,hjust=0,vjust=-0.6) -> p
      }
      # down
      if(nrow(Down)>0)
      { 
        Down %>% dplyr::arrange(-.[[colnames(Up)[2]]]) %>% dplyr::slice(1:GeneNameN) -> GeneNameUpp
        Down %>% dplyr::arrange(.[[colnames(Up)[3]]]) %>% dplyr::slice(1:GeneNameN) -> GeneNameUpfc
        rbind(GeneNameUpp,GeneNameUpfc) %>% unique -> GeneNameDown
        p + geom_text(aes(x=GeneNameDown[,3],y=GeneNameDown[,2]),
                      label=GeneNameDown[,4],size=GeneNameSize,hjust=1,vjust=-0.6) -> p
      }
    }else
      {
        GeneNameList %>% as.character %>% unique %>% data.frame(Symbol=.) %>% 
          dplyr::inner_join(dat, by="Symbol") %>% dplyr::select(c(2,3,4,1)) -> GeneName0
        p + geom_point(aes(x=GeneName0[,3],y=GeneName0[,2]),colour="darkorange") -> p
        GeneName0 %>% dplyr::filter(.[[colnames(GeneName0)[3]]] > 0 ) -> GeneNameUp
        GeneName0 %>% dplyr::filter(.[[colnames(GeneName0)[3]]] < 0 ) -> GeneNameDown
        # Up
        if(nrow(GeneNameUp) > 0)
        {
          p + geom_text(aes(x=GeneNameUp[,3],y=GeneNameUp[,2]),label=GeneNameUp[,4],
                        size=GeneNameSize,fontface=2,hjust=-0.2) -> p
        }
        # Down
        if(nrow(GeneNameDown) > 0)
        {
          p + geom_text(aes(x=GeneNameDown[,3],y=GeneNameDown[,2]),label=GeneNameDown[,4],
                        size=GeneNameSize,fontface=2,hjust=1) -> p
        }
        # # up
        # GeneName0[order(-GeneName0[,2]),] %>% dplyr::slice(1:GeneNameN) -> GeneNameUp
        # GeneNameUp %>% slice(which(GeneNameUp[,3] > 1)) -> GeneNameUp
        # p + aes(x=GeneNameUp[,3], y=GeneNameUp[,2]) -> p
        # p + geom_point(aes(x=GeneNameUp[,3], y=GeneNameUp[,2]), colour="red") -> p
        # p + geom_text(aes(x=GeneNameUp[,3], y=GeneNameUp[,2]), label=GeneNameUp[,4], size=GeneNameSize, fontface=2, hjust=-0.2) -> p
        # # down
        # GeneName0[order(-GeneName0[,2]),] %>% dplyr::slice(1:GeneNameN)-> GeneNameDown
        # GeneNameDown %>% dplyr::slice(which(GeneNameDown[,3] < -1)) -> GeneNameDown
        # p + geom_point(aes(x=GeneNameDown[,3], y=GeneNameDown[,2]), colour="green4") -> p
        # p + geom_text(aes(x=GeneNameDown[,3], y=GeneNameDown[,2]), label=GeneNameDown[,4], size=GeneNameSize, fontface=2, hjust=1) -> p
      }
  }
  p + theme_bw() -> p
  p + theme(plot.title=element_text(size=10, face="bold", hjust=0.5), legend.position="none") -> p
  p + geom_hline(yintercept=-log10(pval), linetype="solid", colour="orange", size=0.1) -> p
  p + geom_vline(xintercept=log2(fc), linetype="solid", colour="orange", size=0.1) -> p
  p + geom_vline(xintercept=-log2(fc), linetype="solid", colour="orange", size=0.1 ) -> p
  p + labs(x="log2ratio(fc)", y="-log10(p)") -> p
  p + ggtitle(paste(title, sep="")) -> p
  p %>% return()
}
