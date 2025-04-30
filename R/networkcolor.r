#' @title Set network node color
#' @description set network node color
#' @param foldchange data.frame Symbol list with fold-change values to display
#' @param colornode character vector of color to set
#' @param networkobj igraph An igraph object.
#' @return no values
#' @examples
#' # not run
#' # acp( mat1 , sif1 )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom gplots colorpanel
#' @importFrom stats prcomp
#' @importFrom scales percent
#' @importFrom graphics hist
#' @noRd
networkcolor <- function( foldchange = NULL, colornode = NULL, networkobj = NULL )
{
  j=1
  colornode -> colorNode
  #
  foldchange %>% setNames(c("Symbol","FC")) -> FoldChange0
  FoldChange0$Symbol %>% as.character %>% is.na %>% which -> sel
  if(length(sel)>0){ FoldChange0[-sel,] -> FoldChange0 }
  FoldChange0$Symbol %>% as.character %>% "=="(.,"NA") %>% which -> sel
  if(length(sel)>0){ FoldChange0[-sel,] -> FoldChange0 }
  FoldChange0$FC %>% is.na %>% which -> sel
  if(length(sel)>0){ FoldChange0[-sel,] -> FoldChange0 }
  FoldChange0$FC %>% as.character %>% "=="(.,"NA") %>% which -> sel
  if(length(sel)>0){ FoldChange0[-sel,] -> FoldChange0 }
  FoldChange0 %>% group_by(.data$Symbol) %>% dplyr::slice_max(.data$FC) %>% data.frame -> FoldChange1
  if( V(networkobj)$Symbol %in% FoldChange1$Symbol %>% which %>% length %>% ">"(.,0) ){FC <- TRUE}else{ FC <- FALSE ; FoldChange1 -> FoldChange2 }
  FALSE -> UP
  ###
  if( FC )
  {
    V(networkobj)$Symbol %>% as.character %>% data.frame("Symbol"=.,stringsAsFactors=F) %>% inner_join(FoldChange1) %>% data.frame -> FoldChange2
    # up
    FoldChange2$FC %>% ">="(1) %>% which -> selup
    if(length(selup) > 0)
    {
      FoldChange2[selup,] %>% dplyr::arrange(.data$FC) %>% data.frame -> FoldChange2Up0
      ifelse( FoldChange2Up0 %>% length %>% ">"(.,50), Breaks0 <- 100, Breaks0 <- 25 )
      FoldChange2Up0$FC %>% hist(breaks=Breaks0, plot = FALSE ) -> hUp0
      FoldChange2Up0$FC %>% hist(breaks=Breaks0, plot = T )
      hUp0$counts -> hUp1
      hUp1
      keycolor = c("gold","orange","red")
      gplots::colorpanel(length(hUp1), low=keycolor[1], mid=keycolor[2], high=keycolor[3]) -> colorFCUp0
      hUp1 %>% "=="(.,0) %>% "!"(.) %>% which -> sel
      if(length(sel) > 0)
      {
        hUp1[sel] -> hUp2
        colorFCUp0[sel] -> colorFCUp1
      }else{
        hUp1 -> hUp2
        colorFCUp0 -> colorFCUp1 }
      hUp3 <- foreach(j=1:length(hUp2), .combine = "c") %do% { sum(hUp2[1:j]) }
      hUp4 <- foreach(j=1:length(hUp3)) %do% { list( rep( colorFCUp1[j], hUp2[j] ), FoldChange2Up0$FC[ (hUp3[j]-(hUp2[j]-1)):hUp3[j] ] ) }
      hUp4 %>% lapply("[[",1)  %>% unlist -> hUp5
      V(networkobj)$Symbol %in% FoldChange2Up0$Symbol %>% which -> sel
      if( length(sel) > 0 ){ colornode %>% replace(sel,hUp5) -> colorNode ; TRUE -> UP }
    }
    # down
    FoldChange2$FC %>% "<="(-1) %>% which -> sel
    if(length(sel)>0)
    {
      FoldChange2[sel,] %>% dplyr::arrange( abs(.data$FC) ) %>% data.frame -> FoldChange2Down0 ; FoldChange2Down0$FC <- abs(FoldChange2Down0$FC)
      ifelse( FoldChange2Down0 %>% length %>% ">"(.,50), Breaks0 <- 100, Breaks0 <- 25 )
      FoldChange2Down0$FC %>% hist(breaks=Breaks0, plot = FALSE ) -> hDown0
      hDown0$counts -> hDown1
      keycolor = c( "greenyellow", "green" , "green4")
      colorpanel(length(hDown1), low=keycolor[1], mid=keycolor[2], high=keycolor[3]) -> colorFCDown0
      hDown1 %>% "=="(.,0) %>% "!"(.) %>% which -> sel
      if( length(sel)>0 )
      {
        hDown1[sel] -> hDown2
        colorFCDown0[sel] -> colorFCDown1
      }else{
        hDown1 -> hDown2
        colorFCDown0 -> colorFCDown1 }
      hDown3 <- foreach(j=1:length(hDown2), .combine = "c") %do% { sum(hDown2[1:j]) }
      hDown4 <- foreach(j=1:length(hDown3)) %do% { list( rep( colorFCDown1[j], hDown2[j] ), FoldChange2Down0$FC[ (hDown3[j]-(hDown2[j]-1)):hDown3[j] ] ) }
      hDown4 %>% lapply("[[",1)  %>% unlist -> hDown5
      V(networkobj)$Symbol %in% FoldChange2Down0$Symbol %>% which -> sel
      if( length(sel) > 0 & !UP ){ colornode %>% replace(sel,hDown5) -> colorNode }
      if( length(sel) > 0 & UP ){ colorNode %>% replace(sel,hDown5) -> colorNode }
    }
  }
  list(colorNode,FoldChange2) %>% return()
}
