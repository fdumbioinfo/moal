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
#' # anovapatternfilter( dat )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select slice filter
#' @importFrom foreach foreach %do%
#' @noRd
anovapatternfilter2fc <- function(
    dat , sif , annot , comp , threshold = c(c(2,6,11),c(2,6,11)+120),
    path = "." , dirname = "DiffLists" )
{
  i=j=1
  thresholdlist[threshold] -> Threshold0
  # colFactor
  comp %>% paste("^p_",.,"$",sep="") %>% paste0(collapse="|") -> grep
  dat %>% colnames %>% grep(grep,.) -> ColFactor
  # col0
  comp %>% paste("^",.,"$",sep="") %>% paste0(collapse="|") -> grep
  sif %>% colnames %>% grep(grep,.) %>% sif[,.] -> Comp0
  Comp0 %>% levels %>% rev %>% as.character %>% combn(2,simplify=F) %>%
    lapply(paste0,collapse="vs") %>% unlist %>% paste("^p_",.,"$",sep="") %>%
    paste0(collapse="|") %>% paste(.,sep="") -> grep
  dat %>% colnames %>% grep(grep,.,value=F) -> Col0
  Comp0 %>% levels %>% length -> NbLevelFactor0
  if(NbLevelFactor0 > 2)
  {
    path %>% file.path(dirname) -> Path ; Path %>% dir.create
    # i=j=1
    j=1
    Col1 <- foreach(i=1:(NbLevelFactor0-1) ) %do%
      {
        j:(j+(NbLevelFactor0-1)-i) -> sel
        j + (NbLevelFactor0-i) -> j
        Col0[sel]
      }
    # i=j=1
    j=1
    Col2 <- foreach( i = 1:length( Col1 ) ) %do%
      {
        foreach(j=1:(length(Col1[[i]])),.combine="c") %do%
          { Col1[[i]] %>% as.character %>% combn( j , simplify = F ) }
      }
    j=1
    ColPattern1 <- foreach(i=1:(length(Col1)-1)) %do%
      { foreach(j=1:length(Col1[[i]]),.combine="c") %do% { Col1[[j]][1] %>% as.character  } }
    Col2 %>% unlist(recursive=F) %>% c(ColPattern1) -> Col3
    # i=1
    #
    foreach(i=1:length(Threshold0)) %do%
    {
      paste("p",Threshold0[[i]][[3]],"_fc",Threshold0[[i]][[4]],"_",as.character(comp),sep="") -> DirNameThr
      Path %>% file.path(DirNameThr) %>% dir.create
      ### Comparisons p-value and fold-change
      # paste("pattern",sep="") -> DirNamePattern
      # Path %>% file.path( DirNameThr , DirNamePattern ) %>% dir.create
      # # j=1
      # #
      # foreach(j=1:length(Col3)) %do%
      #   {
      #     k=1
      #     Patterns0 <- foreach(k=1:length(Col3[[j]])) %do% { c(">","< - ") }
      #     Patterns0 %>% expand.grid -> Patterns1
      #     l=1
      #     foreach(l=1:nrow(as.data.frame(Patterns1))) %do%
      #       {
      #         m=1
      #         Sel0 <- foreach(m=1:ncol(as.data.frame(Patterns1))) %do%
      #           {
      #             paste0(
      #               ".[[ as.numeric( Col3[[j]][m] ) ]] ", "< Threshold0[[i]][[1]] & ",
      #               ".[[ as.numeric( Col3[[j]][m] ) + 1 ]] ", as.character( as.data.frame( Patterns1 )[l,m] ),
      #               " Threshold0[[i]][[2]]", "") -> Text
      #             dat %>% filter( eval( parse( text = Text ) ) ) %>% dplyr::select( c(1) ) %>% unlist %>% as.character
      #           }
      #         dat %>% dplyr::select(c( 1,c(as.numeric(Col3[[j]])+1,as.numeric(Col3[[j]])) %>% sort)) %>% cbind(annot[,-1]) -> r1
      #         Sel0 %>% Reduce(intersect,.) %>% data.frame( "rowID"=. , stringsAsFactors=F ) %>% inner_join(r1,by="rowID" ) %>% arrange(.[[2]]) -> r2
      #         as.character(unlist(Patterns1[l,])) %>% unlist %>% as.character %>% paste0(collapse="") %>% gsub(">","u",.) %>%
      #           gsub("< - ","d",.) %>% paste("List_",DirNameThr,"_",.,"_",sep="") %>%
      #           paste(colnames(dat)[as.numeric(Col3[[j]])] %>% sub("p_(.*)","\\1",.) %>% strsplit("vs") %>% unlist %>% unique(fromLast=T) %>% paste0(collapse = "x"),sep="") %>%
      #           paste( "_" , nrow( r2 ) , ".tsv" , sep = "" ) -> FileName
      #         Path %>% file.path(DirNameThr,DirNamePattern,FileName) %>% output(r2,.)
      #       }
      #   }
      ### Comparisons fold-change
      paste("patternFC",sep="") -> DirNamePatternFC
      Path %>% file.path( DirNameThr , DirNamePatternFC ) %>% dir.create
      j=1
      foreach(j=1:length(Col3)) %do%
        {
          k=1
          Patterns0 <- foreach(k=1:length(Col3[[j]])) %do% { c(">","< - ") }
          Patterns0 %>% expand.grid -> Patterns1
          l=1
          foreach(l=1:nrow(as.data.frame(Patterns1))) %do%
            {
              m=1
              Sel0 <- foreach(m=1:ncol(as.data.frame(Patterns1))) %do%
                {
                  paste0(
                    # pval anova
                    ".[[ as.numeric( ColFactor ) ]] ", "< Threshold0[[i]][[1]] & ",
                    # FC
                    ".[[ as.numeric( Col3[[j]][m] ) + 1 ]] ", as.character( as.data.frame( Patterns1 )[l,m] ),
                    " Threshold0[[i]][[2]]", "") -> Text
                  dat %>% dplyr::filter(eval(parse(text=Text))) %>% dplyr::select(c(1)) %>% unlist %>% as.character
                }
              dat %>% dplyr::select(c(1,as.numeric(ColFactor),c(as.numeric(Col3[[j]])+1,as.numeric(Col3[[j]])) %>% sort)) %>% cbind( annot[,-1] ) -> r1
              Sel0 %>% Reduce(intersect,.) %>% data.frame("rowID"=.,stringsAsFactors=F) %>% inner_join(r1,by="rowID") %>% arrange(.[[2]]) -> r2
              as.character(unlist(Patterns1[l,])) %>% unlist %>% as.character %>% paste0(collapse="") %>% gsub(">","u",.) %>% gsub("< - ","d",.) %>%
                paste( "List_" , DirNameThr , "_" , . , "_"  ,sep = "" ) %>%
                paste(colnames(dat)[as.numeric(Col3[[j]])] %>% sub("p_(.*)","\\1",.) %>% strsplit("vs") %>% unlist %>% unique(fromLast=T) %>% paste0(collapse="x"),sep="") %>%
                paste("_",nrow(r2),".tsv",sep="") -> FileName
              Path %>% file.path(DirNameThr,DirNamePatternFC,FileName) %>% output(r2,.)
            }
        }
      ### sort lists in pattern directory
      # remove 1 comp pattern and 0 list
      Path %>% file.path(DirNameThr) %>% list.files(recursive=T,full.names=T) -> l0
      l0 %>% basename %>% grep("_u_|_d_",.) -> sel
      l0 %>% basename %>% sub(".*_(.*).tsv","\\1",.) %>% as.numeric %>% "=="(.,0) %>% which %>% "c"(.,sel) %>% unique -> sel
      if(length(sel) > 0){l0[sel] %>% file.remove }
      # remove vs Ctrl patterns
      Path %>% file.path(DirNameThr) %>% list.files(recursive=T,full.names=T) -> l0
      Comp0 %>% levels %>% "["(1) %>% paste(".*",.,"_(.*)\\.tsv$",sep="") -> GrepRemove
      GrepRemove
      l0 %>% basename %>% grep(GrepRemove,.) -> sel
      if(length(sel) > 0){l0[sel] %>% file.remove }
      #
      Path %>% file.path(DirNameThr) %>% list.files(recursive=T,full.names=T) -> l0
      #
      if( length(l0) > 0 )
      {
        # pattern
        # Path %>% file.path(DirNameThr,"pattern") %>% list.files( recursive = T , full.names = T ) -> l0
        # if(length(l0) > 0)
        # {
        #   l0 %>% basename %>% gsub("_FC","",.) %>% strsplit("_") %>% lapply("[[",5) %>% unlist -> ll0
        #   ll0 %>% nchar -> ll1
        #   j=1
        #   foreach(j=1:length(table(ll1))) %do%
        #     {
        #       Path %>% file.path(DirNameThr,"pattern",unique(ll1)[j] ) %>% dir.create
        #       ll1 %>% "=="(.,unique(ll1)[j]) %>% which %>% l0[.] -> ll2
        #       ll2 %>% file.copy(Path %>% file.path( DirNameThr , "pattern" , unique(ll1)[j]) , overwrite = T)
        #       ll2 %>% file.remove
        #     }
        # }
        # patternFC
        Path %>% file.path( DirNameThr , "patternFC") %>% list.files( recursive = T , full.names = T ) -> l0
        l0 %>% basename %>% gsub("_FC" , "" , . ) %>% strsplit("_") %>% lapply("[[" , 5 ) %>% unlist -> ll0
        ll0 %>% nchar -> ll1
        if(length(l0)>0)
        {
          l0 %>% basename %>% gsub("_FC","",.) %>% strsplit("_") %>% lapply("[[",5) %>% unlist -> ll0
          ll0 %>% nchar -> ll1
          j=1
          foreach(j=1:length(table(ll1))) %do%
            {
              Path %>% file.path(DirNameThr,"patternFC",unique(ll1)[j] ) %>% dir.create
              ll1 %>% "=="(.,unique(ll1)[j]) %>% which %>% l0[.] -> ll2
              ll2 %>% file.copy(Path %>% file.path( DirNameThr , "patternFC" , unique(ll1)[j]) , overwrite = T)
              ll2 %>% file.remove
            }
        }
      }

    }
  }
}