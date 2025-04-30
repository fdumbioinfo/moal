#' @title anova statistics
#' @param dat data.frame from anova function
#' @param path character
#' @param dirname character
#' @details
#' The function will create a directory including 3 plots of anova statistics computing to obtain pvalues:
#' SumSq2 , MeanSumSq2 and Fratio barplot 
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select contains
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom foreach %dopar% %do%
#' @noRd
anovastats <- function( dat , path = "." , filename = "anovastatplot" )
{
  # file.path(path,dirname) -> Path
  # if(!file.exists(Path)){Path %>% dir.create}
  # output
  dat[,-1] -> Dat1
  # paste("anovastats_data_",ncol(Dat1),"_",nrow(Dat1),".tsv",sep="") -> FileName
  # Dat1 %>% data.frame( rowid = dat[,1] , . , stringsAsFactors = FALSE ) %>% output( file.path( Path , FileName ) )
  # Factor names
  Dat1 %>% colnames %>% sub(".*_(.*)","\\1",.) %>% unique %>% as.factor -> Factor
  #
  # Sum of square
  #
  # Dat1 %>% data.frame %>% dplyr::select(contains("Sumsq")) %>% apply(2 , mean ) %>%
  #   "/"( . , sum(.) ) %>% "*"( . , 100 ) %>% data.frame( SumSq2 = . , Factor ) -> Dat2
  # Dat2 %>%
  #   ggplot( aes( x = "" , y = .data$SumSq2 , fill = Factor ) ) +
  #   geom_bar(  stat = "identity" ) +
  #   scale_fill_manual(values = palette0 ) +
  #   ggtitle( paste( "Average sum of squares" , sep = "" ) ) +
  #   coord_polar( theta = "y") +
  #   theme_bw() +
  #   theme( axis.text.x = element_text( face = "plain", color = "black", size = 7, angle = 0) ) +
  #   xlab( ""  ) + ylab("")
  # paste("Average_Sum_of_squares_",nrow(dat),".pdf" , sep = "" ) -> FileName
  # ggsave( file.path(Path , FileName ) )
  #
  # Mean of square
  #
  # Dat1 %>% data.frame %>% dplyr::select( contains( "Meansq" ) ) %>% apply( 2 , mean ) %>%
  #   "/"( . , sum(.) ) %>% "*"( . , 100 ) %>% data.frame( Meansq = . , Factor ) -> Dat2
  # Dat2 %>%
  #   ggplot( aes( x = "" , y = .data$Meansq , fill = Factor ) ) +
  #   geom_bar(  stat = "identity" ) +
  #   scale_fill_manual(values = palette0 ) +
  #   ggtitle( paste( "Average mean of squares" , sep = "" ) ) +
  #   coord_polar( theta = "y") +
  #   theme_bw() +
  #   theme( axis.text.x = element_text( face = "plain", color = "black", size = 7, angle = 0) ) +
  #   xlab( ""  ) + ylab("")
  # paste("Average_mean_of_squares_",nrow(dat),".pdf" , sep = "" ) -> FileName
  # ggsave( file.path(Path , FileName ) )
  #
  # F ratio
  #
  Dat1 %>% data.frame %>% dplyr::select(contains("Fratio")) %>% apply(2,mean) %>% data.frame(Fratio=.,Factor) -> Dat2
  Dat2 %>%
    ggplot(aes(x=Factor,y=.data$Fratio,fill=Factor)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=palette0) +
    theme_bw() +
    ggtitle(paste("Average F ratio ",sep="")) +
    theme(axis.text.x=element_text(face="plain",color="black",size=7,angle=90))
  paste("Average_Fratio_",filename,"_",nrow(dat),".pdf",sep="") -> FileName0
  ggsave(file.path(path,FileName0))
}