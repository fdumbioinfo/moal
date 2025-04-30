#' Function to filter anova results
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
anovafilter <- function(
  dat , sif , annot,
  comp , threshold = 1:25,
  path = "." , dirname = "DiffLists" )
{
  path %>% file.path(dirname) -> Path
  if(!file.exists(Path)){Path %>% dir.create}
  thresholdlist[threshold] -> Threshold0
  ### colFactor
  comp %>% paste( "^p_" , . , "$" , sep = ""  ) %>% paste0( collapse = "|" ) -> grep
  dat %>% colnames %>% grep( grep , . ) -> ColFactor
  ### col0
  comp %>% paste( "^" , . , "$" , sep = "" ) %>% paste0( collapse = "|" ) -> grep
  sif %>% colnames %>% grep( grep , . ) %>% sif[ , . ] -> Comp0
  Comp0 %>% levels %>% rev %>% as.character %>% combn( 2 , simplify = F ) %>%
    lapply( paste0 , collapse = "vs" ) %>% unlist %>% paste( "^p_" , . , "$",  sep = "" ) %>%
    paste0( collapse = "|" ) %>% paste( . , sep = "" ) -> grep
  dat %>% colnames %>% grep( grep , . , value = F ) -> Col0
  Comp0 %>% levels %>% length -> NbLevelFactor0
  Threshold0 %>% sapply( "[[" , c(1) ) %>% unique -> Threshold1
  Threshold0 %>% sapply( "[[" , c(3) ) %>% unique -> Threshold2
  i = 1
  foreach( i=1:length(Threshold1), .packages = c("magrittr","moal","dplyr") ) %dopar%
    {
      paste("p",Threshold2[i],"_fc1_",as.character(comp),sep="") -> DirNameThr
      Path %>% file.path( DirNameThr ) %>% dir.create
      dat %>% filter( .[[ as.numeric( ColFactor ) ]] < Threshold1[i] ) %>%
        select(c(1,as.numeric(ColFactor),as.numeric(Col0), as.numeric(Col0+1))) %>% inner_join( annot , by = "rowID" ) -> r1
      paste("List_p",Threshold2[i],"_fc1_",colnames(dat)[ColFactor] %>% sub("p_","",.),"_",nrow(r1),".tsv",sep="") -> FileName
      Path %>% file.path( DirNameThr , FileName ) %>% output( r1 , . )
    }
}