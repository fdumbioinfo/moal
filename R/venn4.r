#' @title Venn diagramm of 4 lists
#' @description To make a Venn diagramm of 4 lists
#' @param list list of 4 character vector or list of two data.frame to compare
#' @param listnames character list names to display on graph
#' @param returnlist logical
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
#' #   letters[1:15] ,
#' #   c( letters[2:5] , letters[8:23] ) ) %>% venn4
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import ggplot2
#' @importFrom gplots venn
#' @importFrom ggforce geom_ellipse
#' @importFrom rlang .data
#' @importFrom dplyr full_join inner_join
#' @importFrom utils combn
#' @importFrom foreach %dopar% %do%
#' @noRd
venn4 <- function(
  list = NULL , listnames = NULL , returnlist = F,
  title = "Venn Diagram" , plot = T,
  export = F, path = "." , dirname = NULL )
{
  i=1
  list %>% lapply(as.character) -> list
  # listnames
  listnames -> ListNames
  if( is.null( listnames ) )
  { c( 'A','B','C','D' ) %>% paste(. , " (", list %>% sapply(length) , ")" , sep ="" ) -> labels }else{
    listnames %>% paste( c("A: ","B: ","C: ", "D: "),.," (",list %>% sapply(length) ,")",sep="") -> labels }
  # names for case of intersection with 0 values
  names(list) <- LETTERS[1:4]
  list %>% gplots::venn( data = . ,show.plot = F) %>% attr(. , "intersections") -> list1
  c( "A","B","C","D",
     "A:B","A:C","A:D","B:C","B:D","C:D",
     "A:B:C","A:B:D","A:C:D","B:C:D","A:B:C:D" ) -> ListNames
  rep("",15) %>% setNames( ListNames ) -> ListTemp
   foreach( i=1:length(list1) ) %do% { names(list1[i]) %>% "=="(.,names( ListTemp ) ) %>% which -> sel ; ListTemp[sel] <- list1[i] }
  ListTemp -> list2
  list1 %>% names %>% match( ListNames ) %>% replace( rep(0,15) , . , list1 %>% sapply( length ) ) -> listSize
  # coordinate for overlap values
  data.frame(
    listSize,
    x = c(
      -7.5 , -3.0 , +3.0 , +7.5, #spec
      -5.0 , -4.0 , +0.0 , +0.0 , +4.0 , +5.0, #comb2
      -2.3 , +2.0 , -2.0 , +2.3, #comb3
      +0.0 ),
    y = c(
      +2.0 , +6.5 , +6.5 , +2.0, #spec
      +4.0 , -2.5 , -5.0 , +4.0 , -2.5 , +4.0, #comb2
      +1.4 , -3.3 , -3.3 , +1.4, #comb3
      -1.0 ) ) -> v1
  # dt for the legend
  data.frame(
    x = c( -3.60 , -1 , 1 , 3.60 ),
    y = c( 0 , 2 , 2 , 0 ),
    a = rep( 8 , 4 ),
    b = rep( 4 , 4 ),
    angle = c( -pi/4 , -pi/4 , pi/4 , pi/4 ),
    labels = labels ) -> df.venn
  # venn display
  ggplot(df.venn) -> p
  p + ggforce::geom_ellipse(
    aes(
      x0 = .data$x, y0 = .data$y,
      a = .data$a , b = .data$b , angle = .data$angle , fill = labels ),
    alpha = .3, linewidth = 1, colour = 'black' ) -> p
  p + coord_fixed() -> p
  p + theme_void() -> p
  p + theme(
    legend.position = 'bottom',
    legend.text = element_text( size = 7 ),
    legend.key.size  = unit( 1.5 , "line") ) -> p
  p + scale_fill_manual( values = palette0[1:4] ) -> p
  p + scale_colour_manual( values = palette0[1:4] , guide = "none") -> p
  p + labs(fill = NULL) -> p
  p + guides( fill = guide_legend( nrow = 2 , byrow = TRUE ) ) -> p
  p + annotate("text", x = v1$x, y = v1$y, label = v1$listSize, size = 8 ) -> p
  p + ggtitle( title ) -> p
  p + theme( plot.title = element_text(size=7) ) -> p
  # export
  if( export )
  {
    if( is.null(dirname) ){ paste("venn4_",listSize[15], sep="") -> DirName }else{
      paste( dirname,"_",listSize[15] , sep = "" ) -> DirName }
    path %>% file.path( DirName ) %>% dir.create
    c( "SpeA","SpeB","SpeC","SpeD",
       "SpeIntAB","SpeIntAC","SpeIntAD","SpeIntBC","SpeIntBD","SpeIntCD",
       "IntABC","IntABD","IntACD","IntBCD","IntABCD") -> ListNamesOutput0
    listSize %>% ">"(.,0) %>% which -> sel
    ListNamesOutput0[sel]-> ListNamesOutput1
    list2[sel]-> list3
    foreach( i=1:length(list3) ) %do%
      {
        paste("List_",ListNamesOutput1[i],"_",length(list1[[i]]),".tsv",sep="") -> FileName0
        path %>% file.path( DirName , FileName0 ) -> FileName1
        data.frame( "rowID" = list1[[i]] %>% as.character ) %>% output(FileName1)
      }
    # plot
    paste( "Venn4_",listSize[15],".pdf" , sep = "" ) -> FileName
    path %>% file.path( DirName , FileName ) -> FileName
    ggsave( FileName , plot = p )
  }
  if( plot ){ p %>% print }
}