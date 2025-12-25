#' @title Venn diagramm for 2 lists
#' @description To make a Venn diagramm of 3 lists
#' @param list 2 vector or data.frame
#' @param listnames character list names
#' @param returnlist logical
#' @param title character title to display on graph
#' @param plot logical to display the plot or not
#' @param export logical export  results in directory
#' @param path character
#' @param dirname character name of the directory created when export = T
#' @details
#' To make a Venn diagramm of 3 lists.
#' For data.frame the first column of each dt will be taken for comparison
#' Export only for data.frame list
#' @return list list of the 7 sectors of the graph
#' @examples
#' # not run
#' # venn2( list( v1, v2 ) )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @importFrom rlang .data
#' @importFrom dplyr full_join inner_join
#' @noRd
venn2 <- function(
    list = NULL , listnames = NULL , returnlist = F,
    title = "Venn Diagram", plot = T , export = F,
    path = "." , dirname = NULL )
{
  i=1
  mode(list[[1]]) <- "character"
  mode(list[[2]]) <- "character"
  # create labels for legend
  if( is.null(listnames) ) { c( 'A','B' ) %>% paste( .," (",list %>% sapply(length),")", sep ="" ) -> labels}else{
      listnames %>% paste( c("A: ","B: "),.," (", list %>% sapply(length) ,")", sep ="" ) -> labels }
  # names for case of intersection with 0 values
  # Int
  list[[1]] %>% "%in%"( list[[2]] ) %>% which %>% "["( list[[1]] , . ) -> IntAB
  # Spe A
  list[[1]] %>% "%in%"( list[[2]] ) %>% "!"(.) %>% which %>% "["( list[[1]] , . ) -> SpeA
  # Spe A
  list[[2]] %>% "%in%"( list[[1]] ) %>% "!"(.) %>% which %>% "["( list[[2]] , . ) -> SpeB
  c("A","B","AB") -> IntNames
  list( SpeA,SpeB,IntAB ) %>% sapply(length) -> NIntall
  data.frame(
    NIntall,
    x = c( -1.2 , +1.2 , +0.0 ),
    y = c( -0.5 , -0.5 , -0.5 ) ) -> v1
  # dt for the legend
  data.frame(
    x = c(  -0.866 , +0.866 ),
    y = c(  -0.500 , -0.500 ),
    labels = labels ) -> df.venn
  # venn display
  ggplot(df.venn) -> p
  p + ggforce::geom_circle(aes( x0 = .data$x, y0 = .data$y, r = 1.5 , fill = labels), alpha = .3, linewidth = 1, colour = 'black' ) -> p
  p + coord_fixed() -> p
  p + theme_void() -> p
  p + theme( legend.position = 'bottom', legend.text = element_text(size = 8) ) -> p
  p + scale_fill_manual( values = palette0[1:2] ) -> p
  p + scale_colour_manual( values = palette0[1:2] , guide = "none" ) -> p
  p + labs(fill = NULL) -> p
  p + annotate("text", x = v1$x, y = v1$y, label = v1$NIntall, size = 8) -> p
  p + ggtitle( title ) -> p
  #
  if( export )
  {
    if( is.null(dirname) ){ paste("venn3_",length(IntAB), sep="") -> DirName }else{
      paste( dirname,"_",length(IntAB) , sep = "" ) -> DirName }
    path %>% file.path( DirName ) %>% dir.create
    list(SpeA,SpeB,IntAB) %>% setNames( c( "SpeA","SpeB","IntAB" ) ) -> List0
    #
    foreach( i=1:length( List0 ) ) %do%
      {
        paste( "List_", names(List0)[i], "_" , length(List0[[i]]),".tsv", sep="" ) -> FileName
        path %>% file.path( DirName , FileName ) -> FileName
        data.frame( "rowID" = List0[[i]] %>% as.character ) %>% output(FileName)
      }
    paste( "Venn2_",length(IntAB),".pdf" , sep = "" ) -> FileName
    path %>% file.path( DirName, FileName ) -> FileName
    ggsave( FileName , plot = p )
  }
  if( plot ){ p %>% print }
}