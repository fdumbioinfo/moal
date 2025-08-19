#' @title Annotation function for Symbol, NCBI or Ensembl IDs 
#' @description Annotation function for Symbol, NCBI or Ensembl IDs 
#' @param symbollist character list of IDs or Symbols
#' @param species character species 'hs' 'mm' 'rn' 'dr' 'ss' (see details for complete list) 
#' @param ortholog logical if TRUE return homo sapiens ortholog for choosen species
#' @param dboutput character database used for output 'ncbi'(default) or 'ebi'
#' @param idtype character database ID accepted: 'SYMBOL'(default), 'GENE', 'ENST', 'ENSG', 'ENSP'
#' @details
#' Use moal:::orthoinfo to see complete species list
#' @return data.frame
#' @examples
#' # not run
#' # annot(EnsemblGeneIDs)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by left_join slice select
#' @importFrom rlang .data
#' @import moalannotgene
#' @import moalannotensg
#' @import moalannotenst
#' @import moalannotensp
#' @export
annot <- function( symbollist, species = NULL, ortholog = F, dboutput = "ncbi",idtype = NULL )
{
  symbollist %>% as.character -> symbollist0
  GeneDb0 <- NULL ; ens0 <- NULL
  #
  # check NA, "", " "
  #
  symbollist0 %>% is.na %>% which -> sel
  if( length(sel)>0 ){ symbollist0 %>% replace( sel, paste("row",sel,sep="") ) %>% as.character -> symbollist0 }
  symbollist0 %>% "=="(.,"") %>% which -> sel
  if( length(sel)>0 ){ symbollist0 %>% replace( sel, paste("row",sel,sep="") ) %>% as.character -> symbollist0 }
  symbollist0 %>% "=="(.," ") %>% which -> sel
  if( length(sel)>0 ){ symbollist0 %>% replace( sel, paste("row",sel,sep="") ) %>% as.character -> symbollist0 }
  #
  # ID check
  #
  if( !is.null(idtype) ){ idtype -> IdType }
  if( is.null(idtype) )
  {
    IdType <- "SYMBOL"
    if( length(symbollist0) > 5 )
    {
      if(symbollist0[1:length(symbollist)] %>% grepl( "^[0-9]+$" , . ) %>% which %>% length %>% ">="(.,(length(symbollist)/2))){ IdType <- "GENE" }
      if(symbollist0[1:length(symbollist)] %>% grepl( "^ENS.*T" , . ) %>% which %>% length %>% ">="(.,(length(symbollist)/2))){ IdType <- "ENST" }
      if(symbollist0[1:length(symbollist)] %>% grepl( "^ENS.*G" , . ) %>% which %>% length %>% ">="(.,(length(symbollist)/2))){ IdType <- "ENSG" }
      if(symbollist0[1:length(symbollist)] %>% grepl( "^ENS.*P" , . ) %>% which %>% length %>% ">="(.,(length(symbollist)/2))){ IdType <- "ENSP" }
      if(symbollist0[1:length(symbollist)] %>% grepl( "^[A-Z][0-9].*[0-9]$" , . ) %>% which %>% length %>% ">="(.,(length(symbollist)/2))){ IdType <- "UNIPROT" }
    }
    if( length(symbollist0) <= 5 )
    {
      if(symbollist0[1:5] %>% grepl( "^[0-9]+$" , . ) %>% which %>% length %>% ">"(.,2)){ IdType <- "GENE" }
      if(symbollist0[1:5] %>% grepl( "^ENS.*T" , . ) %>% which %>% length %>% ">"(.,2)){ IdType <- "ENST" }
      if(symbollist0[1:5] %>% grepl( "^ENS.*G" , . ) %>% which %>% length %>% ">"(.,2)){ IdType <- "ENSG" }
      if(symbollist0[1:5] %>% grepl( "^ENS.*P" , . ) %>% which %>% length %>% ">"(.,2)){ IdType <- "ENSP" }
      if(symbollist0[1:5] %>% grepl( "^[A-Z][0-9].*[0-9]$" , . ) %>% which %>% length %>% ">="(.,2)){ IdType <- "UNIPROT" }
    }
  }
  if(ortholog){ IdType <- "ORTHO"}
  # species check
  ifelse( is.null(species), orthoinfo[[6]] -> Species0, orthoinfo %>% sapply("[[",1) %>% grep(species,.) %>% "[["(orthoinfo,.) -> Species0 )
  #
  # NCBI GENE
  #
  if( IdType == "GENE" )
  {
    paste("moalannotgene::",paste0("genedb", Species0[1])," -> GeneDb0",sep="") -> text
    eval(expr = parse(text = text ) )
    GeneDb0 -> GeneDb1
    symbollist0 %>% data.frame("GeneID"=.) %>% left_join(GeneDb1) %>% dplyr::select(-.data$Syn) -> GeneDb1
    GeneDb1 %>% group_by(.data$GeneID) %>% dplyr::slice(1) %>% data.frame -> GeneDb2
    symbollist0 %>% data.frame( "GeneID" =  . , "InputID" = paste("input",1:length(symbollist0),.,sep = "" ) ) -> symbollist1
    GeneDb2 %>% inner_join(symbollist1)  -> GeneDb3
    symbollist0 %>% match(GeneDb3$GeneID) %>% GeneDb3[.,] -> Annot
  }
  #
  # Ensembl gene
  #
  if( IdType == "ENSG" )
  {
    paste("moalannotensg::",paste0("ensg", Species0[1])," -> ens0",sep="") -> text
    eval(expr = parse(text = text ) )
    symbollist0 %>% data.frame("ENSGID"=.) %>% left_join(ens0) -> ens1
    symbollist0 %>% data.frame( "ENSGID" =  . , "InputID" = paste("input",1:length(symbollist0),.,sep = "" ) ) -> symbollist1
    ens1 %>% inner_join(symbollist1)  -> ens2
    symbollist0 %>% match(ens2$ENSGID) %>% ens2[.,] -> Annot
    Annot$Symbol %>% grep("^$",.) -> sell
    if(length(sell)>0){ Annot$ENSGID[sell] -> Annot$Symbol[sell] }
  }
  #
  # Ensembl transcript
  #
  if( IdType == "ENST" )
  {
    paste("moalannotenst::",paste0("enst", Species0[1])," -> ens0",sep="") -> text
    eval(expr = parse(text = text ) )
    symbollist0 %>% data.frame("ENSTID"=.) %>% left_join(ens0) -> ens1
    symbollist0 %>% data.frame( "ENSTID" =  . , "InputID" = paste("input",1:length(symbollist0),.,sep = "" ) ) -> symbollist1
    ens1 %>% inner_join(symbollist1)  -> ens2
    symbollist0 %>% match(ens2$ENSTID) %>% ens2[.,] -> Annot
    Annot$Symbol %>% grep("^$",.) -> sell
    if(length(sell)>0){ Annot$ENSGID[sell] -> Annot$Symbol[sell] }
  }
  #
  # Ensembl Protein
  #
  if( IdType == "ENSP" )
  {
    paste("moalannotensp::",paste0("ensp", Species0[1])," -> ens0",sep="") -> text
    eval(expr = parse(text = text ) )
    symbollist0 %>% data.frame("ENSPID"=.) %>% left_join(ens0) -> ens1
    symbollist0 %>% data.frame( "ENSPID" =  . , "InputID" = paste("input",1:length(symbollist0),.,sep = "" ) ) -> symbollist1
    ens1 %>% inner_join(symbollist1)  -> ens2
    symbollist0 %>% match(ens2$ENSPID) %>% ens2[.,] -> Annot
    Annot$Symbol %>% grep("^$",.) -> sell
    if(length(sell)>0){ Annot$ENSGID[sell] -> Annot$Symbol[sell] }
  }
  #
  # SYMBOL ID
  #
  if( IdType == "SYMBOL" )
  {
    if( dboutput == "ncbi" )
    {
      paste("moalannotgene::",paste0("genedb", Species0[1])," -> GeneDb0",sep="") -> text
      eval(expr = parse(text = text ) )
      GeneDb0 -> GeneDb1
      GeneDb1 %>% group_by(.data$Symbol) %>% dplyr::slice(1) %>% data.frame -> GeneDb1Symb0
      symbollist0 %>% as.character %>% data.frame("Symbol"=.) %>% left_join(GeneDb1Symb0) %>% data.frame -> GeneDb2
      GeneDb2 %>% dplyr::select(-.data$Syn) -> Annot
      # if not 100% match check synonyms
      MATCHALL <- F
      ifelse( GeneDb2$GeneID %>% is.na %>% "!"(.) %>% which %>% length %>% "=="(nrow(GeneDb2)) , MATCHALL <- T, MATCHALL <- F  )
      if( !MATCHALL )
      {
        # check synonyms
        ( ( GeneDb2$GeneID %>% is.na ) & ( GeneDb2$Symbol %>% grepl("^row",.) %>% "!"(.) ) ) %>% which -> sel
        SYN <- F
        ifelse( sel %>% length %>% ">"(.,0), SYN <- T, SYN <- F )
        if( SYN )
        {
          GeneDb1 %>% dplyr::group_by(.data$Syn) %>% dplyr::slice(1) -> GeneDb1Syn0
          GeneDb2[sel,] -> GeneDb3
          GeneDb3$Symbol %>% as.character %>% "%in%"(GeneDb1Syn0$Syn) %>% which -> selsel
          if( selsel %>% length %>% ">"(.,0) )
          {
            GeneDb3[selsel,] -> GeneDb4
            GeneDb4$Symbol %>% as.character %>% data.frame("Syn"=.) %>% left_join(GeneDb1Syn0) %>% dplyr::select(-.data$Syn) -> GeneDb2Syn1
            GeneDb2Syn1[,c(2,1,3:ncol(GeneDb2Syn1))] -> GeneDb2Syn1   
            GeneDb3[selsel,] <- GeneDb2Syn1
            GeneDb2[sel,] <- GeneDb3
            GeneDb2 %>% dplyr::select(-.data$Syn) -> Annot
          }
        }
      }
    }
    #
    if( dboutput == "ebi" )
    {
      paste("moalannotensg::",paste0("ensg", Species0[1])," -> GeneDb0",sep="") -> text
      eval(expr = parse(text = text ) )
      GeneDb0 -> GeneDb1
      GeneDb1 %>% group_by(.data$Symbol) %>% dplyr::slice(1) %>% data.frame -> GeneDb1Symb0
      symbollist0 %>% as.character %>% data.frame("Symbol"=.) %>% left_join(GeneDb1Symb0) %>% data.frame -> GeneDb2
      GeneDb2 %>% dplyr::select(-.data$Syn) -> Annot
      # if not 100% match check synonyms
      MATCHALL <- F
      ifelse( GeneDb2$ENSGID %>% is.na %>% "!"(.) %>% which %>% length %>% "=="(nrow(GeneDb2)) , MATCHALL <- T, MATCHALL <- F  )
      if( !MATCHALL )
      {
        # check synonyms
        ( ( GeneDb2$ENSGID %>% is.na ) & ( GeneDb2$Symbol %>% grepl("^row",.) %>% "!"(.) ) ) %>% which -> sel
        SYN <- F
        ifelse( sel %>% length %>% ">"(.,0), SYN <- T, SYN <- F )
        if( SYN )
        {
          GeneDb1 %>% dplyr::group_by(.data$Syn) %>% dplyr::slice(1) -> GeneDb1Syn0
          GeneDb2[sel,] -> GeneDb3
          GeneDb3$Symbol %>% as.character %>% "%in%"(GeneDb1Syn0$Syn) %>% which -> selsel
          if( selsel %>% length %>% ">"(.,0) )
          {
            GeneDb3[selsel,] -> GeneDb4
            GeneDb4$Symbol %>% as.character %>% data.frame("Syn"=.) %>% left_join(GeneDb1Syn0) %>% dplyr::select(-.data$Syn) -> GeneDb2Syn1
            GeneDb2Syn1[,c(2,1,3:ncol(GeneDb2Syn1))] -> GeneDb2Syn1   
            GeneDb3[selsel,] <- GeneDb2Syn1
            GeneDb2[sel,] <- GeneDb3
            GeneDb2 %>% dplyr::select(-.data$Syn) -> Annot
          }
        }
      }
    }
  }
  #
  # Ortholog
  #
  if( IdType == "ORTHO" )
  {
    moalannotgene::orthogenedb -> OrthoGeneDb0
    OrthoGeneDb0$Other_tax_id %>% grep(Species0[2],.) -> sel
    OrthoGeneDb0 %>% dplyr::slice( sel ) -> OrthoGeneDb1
    paste("moalannotgene::",paste0("genedb", Species0[1])," -> GeneDb0",sep="") -> text
    eval(expr = parse(text = text ) )
    GeneDb0 -> GeneDb1
    symbollist0 %>% data.frame("Symbol"=.) %>% left_join(GeneDb1) -> GeneDb2Symb
    symbollist0 %>% data.frame("Syn"=.) %>% left_join(GeneDb1) -> GeneDb2Syn
    rbind(GeneDb2Symb,GeneDb2Syn) %>% data.frame -> GeneDb2
    GeneDb2 %>% group_by(.data$Symbol) %>% dplyr::slice(1) %>% 
      data.frame %>% dplyr::select(c(3,2,3:ncol(GeneDb2))) %>% data.frame -> GeneDb3
    GeneDb3$GeneID %>% as.character %>% data.frame( "Other_GeneID" = . ) %>% inner_join(OrthoGeneDb1) %>%
      dplyr::select(.data$GeneID,.data$Other_GeneID) %>% setNames(c("GeneID","OtherGeneID")) -> Annot0
    moalannotgene::genedbhs -> GeneDb1
    Annot0 %>% left_join(GeneDb1) %>% dplyr::select(-.data$Syn) -> GeneDb2
    GeneDb2 %>% group_by(.data$GeneID) %>% dplyr::slice(1) %>% data.frame -> Annot
  }
  Annot %>% return()
}
