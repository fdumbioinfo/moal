#' @title Cluster analysis
#' @description make a cluster analysis on row and cut tree for a range a value
#' @param dat data.frame data
#' @param factor factor for display groups on heatmap
#' @param Symbol character list of Symbols to display on heatmap
#' @param nc max cluster number to cut
#' @param path character relative path to create results directory
#' @param dirname character for output directory name
#' @details
#' make ascendent hierachical clustering analysis on row
#' @return results directory
#' @examples
#' # not run
#' # data(mat1)
#' # data(sif1)
#' # hcluster( dat = mat1 factor = sif1$F1)
#' @importFrom magrittr %>%
#' @importFrom dplyr slice
#' @importFrom graphics abline text
#' @importFrom grDevices pdf graphics.off
#' @importFrom dendextend cutree set get_leaves_branches_col
#' @importFrom foreach %do%
#' @noRd
hcluster <- function( dat , factor = NULL , Symbol = NULL , nc = c(2,3,6,12),
                      path = "." , dirname = "cluster_analysis" )
{
  dat[,-1] -> Dat1
  dirname -> DirName
  path %>% file.path(DirName) -> Path
  Path %>% dir.create
  ### create heatmap with colored cluster
  i=j=1
  if( nc %>% "<"(.,2) %>% any ){ nc %>% "<="(.,2) %>% "!"(.) %>% which %>% nc[.] -> nc }
  foreach(i=nc) %do%
  {
    paste( "cl" , i , sep = ""  ) -> DirName
    Path %>% file.path( DirName ) %>% dir.create
    # heatmap cutree both
    paste("Heatmap_cuttree_",i,"_",dim(Dat1)[2],"_",dim(Dat1)[1],".pdf",sep="") %>% 
      file.path(Path,DirName,.) -> FileName
    pdf(FileName)
    Dat1 %>% heatmap(factor=factor,labCol=colnames(Dat1),cexCol=0.45,labRow=Symbol,k=i)
    graphics.off()
    # # heatmap cutree row
    # paste("Heatmap_cuttree_row_",i,"_",dim(Dat1)[2],"_",dim(Dat1)[1],".pdf",sep="") %>% 
    #   file.path(Path,DirName,.) -> FileName
    # pdf(FileName)
    # Dat1 %>% heatmap(factor=factor,dendrogram = "row",labCol=colnames(Dat1),cexCol=0.45,labRow=Symbol,k=i)
    # graphics.off()
  }
  Dat1 %>% t %>% hc( doplot = F ) -> dend0
  Dend1 <- foreach( i=nc ) %do% { dend0 %>% dendextend::set( "branches_k_color" , value = palette0[1:i] , k=i ) %>% list(.,i) }
  foreach(i=1:length( Dend1 ) ) %do%
  {
    foreach( j = 1:Dend1[[i]][[2]] ) %do%
    {
      Dend1[[i]][[1]] %>% unlist -> Dend2
      Dend1[[i]][[1]] %>% dendextend::get_leaves_branches_col(.) %>%
          "=="( . , palette0[j] ) %>% which %>% Dend2[.] %>% slice( dat , . ) -> cl0
      dirname %>% strsplit("_") %>% unlist %>% "["(c(2:4)) %>% paste0(collapse = "_") %>%
        paste("List_",.,"_cl",palette0[j],"_",dim(cl0)[1],".tsv",sep="") %>%
        file.path( Path , paste("cl",Dend1[[i]][2],sep ="") , . ) -> FileName
      cl0 %>% output( FileName )
    }
  }
}