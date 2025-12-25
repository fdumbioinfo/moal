#' @title Venn diagramm for 3 lists
#' @description To make a Venn diagramm of 3 lists
#' @param list list of 3 character vector or list of two data.frame to compare
#' @param listnames character list names to display on graph
#' @param title character title to display on graph
#' @param plot logical to display the plot or not
#' @param export logical export list in file
#' @param path character
#' @param dirname character name of the directory created when export = T
#' @details To make a Venn diagramm of 3 lists
#' @return list list of the 7 sectors of the graph
#' @examples
#' # library(magrittr)
#' # list(
#' #   c(letters[6:20] , letters[25] ) ,
#' # letters[1:15] ,
#' #  c( letters[2:5] , letters[8:23] ) ) %>% venn3
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @importFrom rlang .data
#' @importFrom dplyr full_join inner_join
#' @noRd
venn3 <- function(
  list = NULL , listnames = NULL , title = "Venn Diagram",
  plot = T , export = F, path = "." , dirname = NULL )
{
  list %>% lapply(as.character) -> list
  listnames -> ListNames
  if( is.null( listnames ) )
  { c( 'A', 'B' , 'C' ) %>% paste( . , " (", list %>% sapply( length ) , ")" , sep ="" ) -> Listnames }else
    { listnames %>% paste(c("A: ","B: ","C: ") , . , " (", list %>% sapply( length ) ,")" , sep ="" ) -> Listnames }
  # names for case of intersection with 0 values
  # Int 1vs2
  list[[1]] %>% "%in%"( list[[2]] ) %>% which %>% "["( list[[1]] , . ) -> IntAB
  # Int 1vs3
  list[[1]] %>% "%in%"( list[[3]] ) %>% which %>% "["( list[[1]] , . ) -> IntAC
  # Int 2vs3
  list[[2]] %>% "%in%"( list[[3]] ) %>% which %>% "["( list[[2]] , . )  -> IntBC
  # Int 1vs2vs3
  IntAB %>% "%in%"( IntAC ) %>% which %>% "["( IntAB , . ) %>% match(IntBC ) %>% "["( IntBC , . ) -> IntABC
  # Spe Int 1vs2
  IntAB %>% "%in%"( IntABC ) %>% "!"(.) %>% which %>% "["( IntAB , . ) -> SpeIntAB
  # Spe Int 1vs3
  IntAC %>% "%in%"( IntABC ) %>% "!"(.) %>% which %>% "["( IntAC , . ) -> SpeIntAC
  # Spe Int 2vs3
  IntBC %>% "%in%"( IntABC ) %>% "!"(.) %>% which %>% "["( IntBC , . ) -> SpeIntBC
  # Spe 1
  list[[1]] %>% "%in%"( c( SpeIntAB, SpeIntAC , IntABC ) ) %>% "!"( . ) %>% which %>% "["( list[[1]] , . ) -> SpeA
  # Spe 2
  list[[2]] %>% "%in%"( c( SpeIntAB, SpeIntBC , IntABC ) ) %>% "!"( . ) %>% which %>% "["( list[[2]] , . ) -> SpeB
  # Spe 3
  list[[3]] %>% "%in%"( c( SpeIntAC, SpeIntBC , IntABC ) ) %>% "!"( . ) %>% which %>% "["( list[[3]] , . ) -> SpeC
  c("A","B","C","AB","AC","BC","ABC") -> IntNames
  list( SpeA,SpeB,SpeC, SpeIntAB,SpeIntAC,SpeIntBC, IntABC ) %>% sapply(length) -> NIntall
  # dt in order to place values in circle
  data.frame(
    NIntall,
    x = c( +0.0 , -1.2 , +1.2 , -0.8 , +0.8 , +0.0 , 0.0 ),
    y = c( +1.2 , -0.6 , -0.6 , +0.5 , +0.5 , -0.8 , 0.0 ) ) -> v1
  # dt for the legend
  data.frame(
    x = c( +0.000 , -0.866 , +0.866 ),
    y = c( +1.000 , -0.500 , -0.500 ),
    labels = Listnames ) -> df.venn
  # display
  ggplot(df.venn) -> p
  p + ggforce::geom_circle(
    aes(
      x0 = .data$x, y0 = .data$y,
      r = 1.5, fill = Listnames),
    alpha = .3, linewidth = 1, colour = 'black' ) -> p
  p + coord_fixed() -> p
  p + theme_void() -> p
  p + theme( legend.position = 'bottom' , legend.text = element_text(size = 6) ) -> p
  p + scale_fill_manual( values = palette0[1:3] ) -> p
  p + scale_colour_manual( values = palette0[1:3] , guide = "none") -> p
  p + labs(fill = NULL) -> p
  p + annotate("text", x = v1$x, y = v1$y, label = v1$NIntall, size = 8 ) -> p
  p + ggtitle(title) -> p
  p + theme(plot.title = element_text(size=7)) -> p
  # export
  if( export )
  {
    if( is.null(dirname) ){ paste("venn3_",length(IntABC), sep="") -> DirName }else{
        paste( dirname,"_",length(IntABC) , sep = "" ) -> DirName }
    path %>% file.path( DirName ) %>% dir.create
    c(SpeIntAB,SpeIntAC,SpeIntBC,IntABC) -> Int2ABC
    c(SpeIntAB,IntABC,SpeIntAC) -> Int2A
    c(SpeIntAB,IntABC,SpeIntBC) -> Int2B
    c(SpeIntAC,IntABC,SpeIntBC) -> Int2C
    list( SpeA,SpeB,SpeC , SpeIntAB,SpeIntAC,SpeIntBC , IntAB,IntAC,IntBC ,
          SpeIntAB,SpeIntAC,SpeIntBC , IntABC,Int2ABC , Int2A,Int2B,Int2C) %>%
      setNames( c( "SpeA","SpeB","SpeC", "SpeIntAB","SpeIntAC","SpeIntBC" , "IntAB","IntAC","IntBC",
                   "SpeIntAB","SpeIntAC","SpeIntBC" , "IntABC","Int2ABC" , "Int2A","Int2B","Int2C" ) ) -> List0
    i = 1
    foreach( i = 1:length( List0 ) ) %do%
    {
      paste("List_",names(List0)[i],"_",length(List0[[i]]),".tsv",sep="") -> FileName
      path %>% file.path( DirName , FileName ) -> FileName
      data.frame( "rowID" = List0[[i]] %>% as.character ) %>% output(FileName)
    }
    paste( "Venn3_",length(IntABC),".pdf" , sep = "" ) -> FileName
    path %>% file.path( DirName , FileName ) -> FileName
    ggsave( FileName , plot = p )
  }
  if( plot ){ p %>% print }
}