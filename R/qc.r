#' @title Quality Controls
#' @description
#' Descriptive analysis applied on column:
#'  - histogram, boxplot, hierarchical clustering and PCA for column
#' @param dat data.frame first column for rowID column  
#' @param sif data.frame sample information file
#' @param dooutputinput logical if TRUE (by default) export input data
#' @param dohisto logical if TRUE (by default) do histogram
#' @param doboxplot logical if TRUE (by default) do boxplot
#' @param dohc logical if TRUE (by default)do hierarchical clustering
#' @param doacp logical if TRUE (by default) do PCA
#' @param breaks numeric break number for histogramm function
#' @param dirname character
#' @param path character
#' @details dat row must be equal to sif row
#' @return directory including analysis pdf plots
#' @examples
#' # not run
#' # qc(dat,sif)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @import gplots
#' @import ggplot2
#' @importFrom foreach foreach %dopar%
#' @importFrom graphics title hist
#' @export
qc <- function(
    dat , sif = NULL, dooutputinput = FALSE,
    dohisto = TRUE, doboxplot = TRUE, dohc = TRUE, doacp = TRUE,
    breaks = 70,
    dirname = NULL , path = "." )
{
  i=1
  dat[,-1] -> Dat1
  if(is.null(dirname)){ paste("QC_",ncol(Dat1),"_",nrow(Dat1),sep="") -> DirName }else{
    paste("QC_",dirname,"_",ncol(Dat1),"_",nrow(Dat1),sep="")-> DirName  }
  path %>% file.path(DirName) -> Path
  Path %>% dir.create
  #
  # output
  #
  if(dooutputinput)
  {
    paste("data_",ncol(Dat1),"_",nrow(Dat1),".tsv",sep="") -> FileName
    Path %>% file.path( FileName ) -> FileName
    data.frame( "rowID" = dat[,1] %>% as.character, Dat1, stringsAsFactors=T ) -> Dat2
    Dat2 %>% output( FileName )
    # metadata
    if(!is.null(sif))
    {
      # for( i in 1:dim(sif)[2] )
      # {
      #   if( ( length( table( as.character(sif[,i]) ) ) < length( as.character( sif[,i] ) ) ) &
      #       length( table( as.character(sif[,i]) ) ) < 28 )
      #   { sif[,i]  <- factor( sif[,i] ) }
      # }
      paste("metadata_",nrow(sif),".tsv",sep="") -> FileName
      Path %>% file.path(FileName) -> FileName
      sif %>% output(FileName)
    }else
      {
        paste("metadata_",ncol(Dat1),".tsv",sep="") -> FileName
        Path %>% file.path( FileName ) -> FileName
        data.frame( SampleID = paste("s",1:ncol(Dat1),sep=""), FileName=colnames(Dat1) ) %>% output(FileName)
      }
  }
  #
  # 1 - histogram
  #
  if(dohisto)
  {
    paste("Histo_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep="") -> FileName
    Path %>% file.path(FileName) -> FileName
    paste("Histogram",ncol(Dat1)," x ",nrow(Dat1)) -> Title
    pdf(FileName)
    Dat1 %>% as.matrix %>% hist(breaks=70, main=Title)
    graphics.off()
  }
  #
  # without sif
  #
  if(is.null(sif))
  {
    # boxplot
    if(doboxplot)
    {
      paste("boxplot_data_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep="") -> FileName
      Path %>% file.path(FileName) -> FileName
      paste("Boxplot data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
      # chunk sample name and mar bottom -> 11
      Dat1 %>% colnames %>% nchar %>% max -> NcharDat0
      Dat1 %>% colnames -> ColNames0
      ColNames0 %>% nchar %>% ">"(30) %>% which -> selC
      ColNames0[selC] %>% strsplit("") %>% lapply("[",1:30) %>% lapply(paste0,collapse="") -> colnames(Dat1)[selC]
      # par
      # 4 -> NcharDat0
      # c(5,30,50) -> Nchar0
      # c(7,10,15) -> Mar0
      # if(NcharDat0 > 50){ Dat1 %>% colnames %>% paste0() }
      # 7 -> MarBot1
      # if(NcharDat0 > Nchar0[1]){NcharDat0 %>% ">="(Nchar0) %>% which %>% max %>% MarBot0[.] -> MarBot1}
      # width with column number
      c(6,10,20,50,150,200) -> Ncol0
      c(7,10,15,20,30,35) -> Width0
      if(ncol(Dat1) < 200){ncol(Dat1) %>% "<="(Ncol0) %>% which %>% min %>% Width0[.] -> Width1}else{
        Width0[length(Width0)] -> Width1 }
      # plot
      pdf(FileName,width=Width1,height=7)
      par(mar=c(11,4,4,2))
      Dat1 %>% graphics::boxplot(outline=F,las=2,cex.axis=0.7)
      title(main=Title)
      graphics.off()
    }
    # hc
    if(dohc)
    {
      paste("hc_data_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
      Path %>% file.path( FileName ) -> FileName
      paste("Hierarchical clustering data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
      # plot
      pdf( FileName )
      Dat1 %>% select(-1) %>% hc( title = Title )
      graphics.off()
    }
    # acp
    if(doacp)
    {
      Dat1 %>% apply(1,sd) %>% "=="(0) %>% which -> selsd
      if( length(selsd) > 0 ){ Dat1[-selsd,] -> t }else{ Dat1 -> t }
      paste("acp_data_",ncol(t),"_",nrow(t),".pdf",sep = "") -> FileName0
      Path %>% file.path(FileName0) -> FileName1
      paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t),sep="") -> Title
      # plot
      t %>% acp(title=Title)
      ggsave(FileName1)
    }
  }
  # with factor
  if(!is.null(sif))
  {
    sif %>% sapply(is.factor) %>% which -> sel
    # boxplot
    if(doboxplot)
    {
      Path %>% file.path("boxplot") %>% dir.create
      foreach(i=1:length(sel)) %do%
        {
          paste("boxplot_data_",colnames(sif)[sel[i] ],"_",ncol(Dat1),"_",nrow(Dat1),".pdf",sep="") -> FileName
          Path %>% file.path("boxplot",FileName) -> FileName
          paste("Boxplot data ",ncol(Dat1)," x ",nrow(Dat1),sep="") -> Title
          # chunk sample name and mar bottom -> 11
          Dat1 %>% colnames %>% nchar %>% max -> NcharDat0
          Dat1 %>% colnames -> ColNames0
          ColNames0 %>% nchar %>% ">"(30) %>% which -> selC
          ColNames0[selC] %>% strsplit("") %>% lapply("[",1:30) %>% lapply(paste0,collapse="") -> colnames(Dat1)[selC]
          # par
          # width with column number
          c(6,10,20,50,150,200) -> Ncol0
          c(7,10,15,20,30,35) -> Width0
          if(ncol(Dat1) < 200){ncol(Dat1) %>% "<="(Ncol0) %>% which %>% min %>% Width0[.] -> Width1}else{
            Width0[length(Width0)] -> Width1 }
          # plot
          pdf(FileName,width=Width1,height=7)
          par(mar=c(11,4,4,2))
          Dat1 %>% boxplot( factor=sif[,sel[i]], outline=F, title=Title, legendtitle=colnames(sif)[sel[i]] )
          graphics.off()
        }
    }
    # hc
    if(dohc)
    {
      Path %>% file.path( "hc") %>% dir.create
      foreach( i=1:length(sel) ) %do%
        {
          paste("hc_data_",colnames(sif)[ sel[i] ],"_",ncol(Dat1),"_",nrow(Dat1),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "hc", FileName ) -> FileName
          paste("Hierarchical clustering data ",ncol(Dat1)," x ",nrow(Dat1),sep = "") -> Title
          #
          pdf( FileName )
          Dat1 %>% hc( factor=sif[,sel[i]], title=Title, legendtitle=colnames( sif )[sel[i]] )
          graphics.off()
        }
    }
    # ACP
    if(doacp)
    {
      # ACP
      Path %>% file.path( "pca") %>% dir.create
      Dat1 %>% apply(1,sd) %>% "=="(0) %>% which -> selsd
      if( length(selsd) > 0 ){ Dat1[-selsd,] -> t }else{ Dat1 -> t }
      # pc 1 2
      foreach( i=1:length(sel) ) %do%
        {
          paste("pca_data_",colnames(sif)[sel[i] ],"_pc12_",ncol(t),"_",nrow(t),".pdf",sep="") -> FileName0
          Path %>% file.path("pca",FileName0) -> FileName1
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t),sep="") -> Title
          #
          t %>% acp(factor=sif[,sel[i]], pc1=1, pc2=2, title=Title, legendtitle=colnames(sif)[sel[i]])
          ggsave(FileName1)
        }
      # pc 1 3
      foreach( i=1:length(sel) ) %do%
        {
          paste("pca_data_",colnames(sif)[ sel[i] ],"_pc13_",ncol(t),"_",nrow(t),".pdf",sep="") -> FileName0
          Path %>% file.path("pca",FileName0) -> FileName1
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t),sep="") -> Title
          #
          t %>% acp(factor=sif[,sel[i]], pc1=1, pc2=3, title=Title, legendtitle=colnames(sif)[sel[i]])
          ggsave(FileName1)
        }
      # pc 1 4
      foreach( i=1:length(sel) ) %do%
        {
          paste("pca_data_",colnames(sif)[ sel[i] ],"_pc14_",ncol(t),"_",nrow(t),".pdf",sep="") -> FileName0
          Path %>% file.path("pca",FileName0) -> FileName1
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t), sep = ""  ) -> Title
          #
          t %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=4, title=Title, legendtitle=colnames(sif)[ sel[i] ] )
          ggsave( FileName1 )
        }
    }
  }
  #
  # with numeric factor
  #
  if(!is.null(sif))
  {
    if( sif %>% sapply(is.numeric) %>% which %>% any)
    {
      sif %>% sapply(is.numeric) %>% which -> sel
      Dat1 %>% apply(1,sd) %>% "=="(0) %>% which -> selsd
      if( length(selsd) > 0 ){ Dat1[-selsd,] -> t }else{ Dat1 -> t }
      # ACP
      if( !dir.exists(Path %>% file.path("pca")) ){ Path %>% file.path("pca") %>% dir.create }
      # pc 1 2
      foreach( i=1:length(sel) ) %do%
        {
          paste("pca_data_",colnames(sif)[sel[i] ],"_pc12_",ncol(t),"_",nrow(t),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "pca", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t),sep="") -> Title
          #
          t %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=2, title=Title, legendtitle=colnames(sif)[sel[i]] )
          ggsave( FileName )
        }
      # pc 1 3
      foreach( i=1:length(sel) ) %do%
        {
          paste("pca_data_",colnames(sif)[ sel[i] ],"_pc13_",ncol(Dat1),"_",nrow(t),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "pca", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t),sep="") -> Title
          #
          t %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=3, title=Title, legendtitle=colnames(sif)[sel[i]] )
          ggsave( FileName )
        }
      # pc 1 4
      foreach( i=1:length(sel) ) %do%
        {
          paste("acp_data_",colnames(sif)[ sel[i] ],"_pc14_",ncol(t),"_",nrow(t),".pdf" , sep = "" ) -> FileName
          Path %>% file.path( "pca", FileName ) -> FileName
          paste("Principal Component Analysis data ",ncol(t)," x ",nrow(t), sep = ""  ) -> Title
          #
          t %>% acp( factor=sif[ , sel[i] ], pc1=1, pc2=4, title=Title, legendtitle=colnames(sif)[ sel[i] ] )
          ggsave( FileName )
        }
    }
  }
}
