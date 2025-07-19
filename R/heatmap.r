#' @title Heatmap
#' @description To make a heatmap
#' @param dat matrix numeric
#' @param factor factor
#' @param method character
#' @param dendrogram character to display 'none', 'row', 'column' or 'both' (by default) dendrograms
#' @param labCol character
#' @param cexCol numeric
#' @param labRow Character
#' @param cexRow numeric
#' @param cexlegend numeric
#' @param k numeric number of clusters to colorize for rows
#' @param keysize numeric
#' @param keycolor character of 3 for low mid high value of the key
#' @param parmar numeric 4 values for margin sizes
#' @param scale numeric standardize row by defaut and column or none accepted
#' @details To make a heatmap from a matrix or a data.frame
#' @return no returned value
#' @examples
#' # not run
#' # library(magrittr)
#' # data(sif1)
#' # data(mat1)
#' # mat1 %>% heatmap(sif1$F3)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom graphics par legend
#' @importFrom gplots heatmap.2 colorpanel
#' @importFrom dendextend set
#' @export
heatmap <- function(
  dat , factor , method = "complete" , dendrogram = "both", k = NULL,
  labCol = "" , cexCol = 0.85 , labRow = "" , cexRow = NULL,
  cexlegend = 0.65 , keysize = 0.9, keycolor = c( "darkgreen", "orange" , "darkred") , parmar = c(5,4,5,6), scale = "row")
{
  par(xpd=T, mar=parmar)
  if(dendrogram == "both")
  {
    dat %>% hc(factor=factor, plot=F , method = method ) -> coldend
    dat %>% t %>% hc(plot=F, method=method) -> rowdend
  }
  if( dendrogram == "none" )
  {
    coldend <- F
    rowdend <- F
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if( dendrogram == "row" )
  {
    coldend <- F
    rowdend <- T
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if(dendrogram == "col")
  {
    coldend <- T
    rowdend <- F
    factor %>% order -> sel
    dat[ , sel] -> dat
    factor[sel] -> factor
    labCol[sel] -> labCol
  }
  if( !is.null(k) ){ rowdend %>% dendextend::set( "branches_k_color", value=palette0[1:k], k=k ) -> rowdend }
  # adapt legend color when not have same group number to display on the heatmap compare to factor to still keep the same color
  factor %>% levels %>% length -> FactorLevelsNumber
  factor %>% as.character %>% table %>% names %>% length -> MatLevelsNumber
  if(FactorLevelsNumber > MatLevelsNumber)
  {
    factor %>% levels -> LevelsLegend0
    palette0[1:length( levels( factor ) ) ] -> ColLegend0
    factor %>% as.character %>% table %>% names %>% "%in%"(factor %>% levels)
    factor %>% levels %>% "%in%"(factor %>% as.character %>% table %>% names) %>% "!"(.) %>% which -> sel
    ColLegend0[-sel] -> ColLegend1
    LevelsLegend0[-sel] -> LevelsLegend1
  }
  #
  # adapt cexCol to matrix size
  # Symbol size
  if(is.null(cexRow))
  {
    c(5,10,20,30,40,50,100,200,300,400,1000) -> NRowheatmap0
    c(1,0.9,0.8,0.7,0.6,0.5,0.35,0.20,0.18,0.15,0.1) -> CexRow0
    if(nrow(dat) < 1000)
      {nrow(dat) %>% "<="(NRowheatmap0) %>% which %>% min %>% CexRow0[.] -> CexRow1}else{ CexRow0[11] -> CexRow1 }
  }else{ cexRow -> CexRow1}
  #
  
  dat %>% as.matrix %>%
    gplots::heatmap.2(
      ColSideColors=factortocolor(factor), labCol=labCol, cexCol=cexCol, srtCol=90, labRow=labRow, cexRow=CexRow1, srtRow = 0,
      trace="none", key=T, keysize=keysize, density.info="none", key.xlab="", key.title=NA,
      dendrogram=dendrogram, Colv=coldend, Rowv=rowdend,
      col=colorpanel(50, low=keycolor[1], mid=keycolor[2], high=keycolor[3]), scale=scale, distfun=dist2 )
  if(FactorLevelsNumber > MatLevelsNumber)
  {
    legend(
      # "topright", inset = c(-0.22,-0.1), legend = levels( factor( factor ) ),
      "topright", inset = c(-0.22,-0.1), legend = LevelsLegend1,
      col = ColLegend1,
      lty = 1, lwd = 5, cex = cexlegend )
  }else{
    legend(
      "topright", inset = c(-0.22,-0.1), legend = levels( factor( factor ) ),
      col = palette0[1:length( levels( factor( factor ) ) ) ],
      lty = 1, lwd = 5, cex = cexlegend )
  }
}
