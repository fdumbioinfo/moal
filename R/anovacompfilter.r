#' @title filter anova results
#' @description
#'  Function to filter anova results after anova analysis.
#' @param dat data.frame of anova results
#' @param sif data.frame sample description including design factors
#' @param annot data.frame description of normalize data variables
#' @param comp character factor name in anova model for pairwise comparisons
#' @param threshold numeric vector from 1 to 15 for p-value and fold-change
#' @param path character path for record results (relative path)
#' @param dirname character results directory name
#' @return
#' directory
#' @examples
#' # not run
#' # anovafilter( dat )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select slice filter
#' @importFrom foreach foreach %do%
#' @noRd
anovacompfilter <- function(
  dat , sif , annot,
  comp , threshold = 1:25,
  path = "." , dirname = "DiffLists" )
{
  path %>% file.path( dirname ) -> Path
  if(!file.exists(Path)){Path %>% dir.create}
  thresholdlist[threshold] -> Threshold0

  # colFactor
  comp %>% paste( "^p_" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
  dat %>% colnames %>% grep( grep , . ) -> ColFactor

  # col0
  comp %>% paste( "^" , . , "$" , sep = "" ) %>% paste0( collapse = "|" ) -> grep
  sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Comp0
  Comp0 %>% levels %>% rev %>% as.character %>% combn( 2 , simplify = F ) %>%
    lapply( paste0 , collapse = "vs" ) %>% unlist %>% paste( "^p_" , . , "$",  sep = "" ) %>%
    paste0( collapse = "|" ) %>% paste( . , sep = "" ) -> grep
  dat %>% colnames %>% grep( grep , . , value = F ) -> Col0
  Comp0 %>% levels %>% length -> NbLevelFactor0

  i = 1
  foreach( i = 1:length( Threshold0 ) , .packages = c("magrittr","dplyr","moal","foreach") ) %dopar%
    {
      paste( "p" , Threshold0[[i]][[3]] , "_fc" , Threshold0[[i]][[4]] , sep = "" ) %>%
        paste( "_" , as.character( comp ) , sep = "" )-> DirNameThr
      Path %>% file.path( DirNameThr ) %>% dir.create
      j = 1
      foreach( j = 1:length( Col0 ) ) %do%
        {
          # up
          dat %>% filter(
            .[[ as.numeric( Col0[j] ) ]] < Threshold0[[i]][[1]] &
            .[[ as.numeric( Col0[j] ) + 1 ]] > Threshold0[[i]][[2]] ) %>%
            select( c( 1 , as.numeric( Col0[j] ) , as.numeric( Col0[j] + 1 ) ) ) %>%
            inner_join( annot , by = "rowID" ) %>% arrange( .[[2]] ) -> r1u

          paste( "List_" , DirNameThr , sep = "" ) %>%
            paste( "_u_" , colnames(dat)[ Col0[j] ] %>% sub( "p_" , "" , . ) , "_" , nrow(r1u) , ".tsv" , sep = "" ) -> FileName
          Path %>% file.path( DirNameThr , FileName ) %>% output( r1u , . )

          # down
          which(
            dat[ , as.numeric( Col0[j] ) ] < Threshold0[[i]][[1]] &
            dat[ , as.numeric( Col0[j] ) + 1 ] < -Threshold0[[i]][[2]] ) -> sel
          dat %>% slice( sel ) %>%
            select( c( 1 , as.numeric( Col0[j] ) , as.numeric( Col0[j]+1 ) ) ) %>%
            inner_join( annot , by = "rowID") %>% arrange( .[[2]] ) -> r1d

          paste( "List_" , DirNameThr , sep = "" ) %>%
            paste( "_d_" , colnames(dat)[ Col0[j] ] %>% sub( "p_" , "" , . ) , "_" , nrow(r1d) , ".tsv" , sep = "" ) -> FileName

          Path %>% file.path( DirNameThr , FileName ) %>% output( r1d , . )

          # up + down
          rbind( r1u , r1d ) -> r2
          paste( "List_" , DirNameThr , sep = "" ) %>%
            paste( "_ud_" , colnames(dat)[ Col0[j] ] %>% sub( "p_" , "" , . ) , "_" , nrow(r2) , ".tsv" , sep = "" ) -> FileName

          Path %>% file.path( DirNameThr , FileName ) %>% output( r2 , . )
        }
    }#rof threshold
}