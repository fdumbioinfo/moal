#' @title replace 0 value by background noise estimate
#' @param dat character vector of mzxml file
#' @param value numeric value to substitute
#' @param factor character factor for remove non random 0 value distribution across factor levels
#' @param ns number 0 number to remove row
#' @param pval number fisher test p-value threshold for 0 distribution accross
#' @param quantv number quantile probability for background noise value (1st percentile by default)
#' @param logbgvalue logical background noise value logged if TRUE
#' @details
#' replace value by column median if group size is one and replace by row group median if group size > 1
#' @return data.frame
#' @examples
#' # not run
#' # replacegroupmed(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach  %do% foreach
#' @importFrom stats median setNames quantile
#' @importFrom graphics hist
#' @importFrom edgeR calcNormFactors DGEList
#' @importFrom stats fisher.test
#' @export
replacena <- function( dat = NULL, value = NULL, factor = NULL, ns = NULL, pval = 0.05, quantv = 0.01, logbgvalue = FALSE )
{
  i=1
  dat %>% dplyr::select(-1) -> Dat0
  dat %>% dplyr::select(1) %>% unlist -> rowID
  Dat0 %>% unlist -> t
  #
  # replace NA by 0
  t %>% is.na %>% which -> selNA1
  t %>% grep("^$",.) -> selNA2
  t %>% "=="(.,"NA") %>% which -> selNA3
  t %>% "<"(.,0) %>% which -> selNA4
  c(selNA1,selNA2,selNA3,selNA4) %>% unique -> selNAall
  if(length(selNAall)>0){ t %>% replace(selNAall,0) %>% matrix(nrow=nrow(Dat0),byrow=F) %>% data.frame -> Dat1 }else{
    t %>% matrix(nrow=nrow(Dat0),byrow=F) %>% data.frame -> Dat1 }
  #
  # if value: replace 0 by input value
  Dat1 %>% unlist -> t
  t %>% "=="(0) %>% which -> sel0
  if(!is.null(value) & length(sel0)>0)
  { 
    t %>% replace(sel0,value) %>% matrix(nrow=nrow(Dat1),byrow=F) %>% 
      data.frame(rowID=dat %>% dplyr::select(1) %>% unlist,.) %>% setNames(dat %>% colnames) -> Dato
  }else{ dat -> Dato }
  #
  # replace 0 by background noise and remove row with 0
  # -----
  # with factor
  # -----
  if(is.null(value) & length(sel0)>0 & !is.null(factor))
  {  
    # find non random 0 row
    Dat1 %>% t %>% data.frame %>% lapply("==",0) %>% lapply(which) %>% lapply(length) %>% unlist -> v0
    # minimum 0 by row to check randam disctribution
    if(is.null(ns))
    {
      factor %>% table %>% lapply("%/%",2) %>% unlist %>% sum -> ns
      if(ns<1){ factor %>% levels %>% length %>% "%/%"(2) -> ns }
    }
    v0 %>% ">"(ns) %>% which -> sel0
    v0 %>% "<="(ns) %>% which -> selkeep
    # chek randam 0 distribution across factor level
    Dat1 %>% dplyr::slice(sel0) -> tt0
    pv <- foreach(i=1:nrow(tt0),.combine="c") %do%
      {
        tt0 %>% dplyr::slice(i) %>% unlist -> v0
        table(factor,v0 == 0) %>% fisher.test %>% unlist %>% "["(1)
      }
    pv %>% "<"(pval) %>% which -> selkeeppval
    if(length(selkeeppval)>0){ sel0[selkeeppval] -> selkeeppval }
    #
    # create data.frame with no 0 and non homogenaous 0 distribition accross sample
    c(selkeep,selkeeppval) %>% unique -> selkeep2
    Dat1 %>% dplyr::slice(selkeep2) %>% data.frame(rowID=rowID[selkeep2],.) -> Dat2
    Dat2 %>% dplyr::select(-1) -> Dat3
    #
    # create background noise matrix
    Dat3 %>% unlist -> t
    t %>% "=="(0) %>% which -> sel00
    # bg value
    if(logbgvalue){ t[-sel00] %>% log2 %>% quantile(probs=quantv,na.rm=T) %>% as.numeric -> q100 }else
      { t[-sel00] %>% quantile(probs=quantv,na.rm=T) %>% as.numeric -> q100 }
    Dat3 %>% edgeR::DGEList() %>% edgeR::calcNormFactors(method = c("TMM")) -> dg0R
    Dat3 %>% t %>% edgeR::DGEList() %>% edgeR::calcNormFactors(method = c("TMMwsp")) -> dg0C
    dg0C$samples$norm.factors %o% dg0R$samples$norm.factors %>% "*"(q100) %>% data.frame -> dtnoise
    #
    # replace 0 by background noise
    Dat3 %>% unlist -> t
    dtnoise %>% unlist -> tnoise
    t %>% unlist %>% "=="(0) %>% which -> sel0
    sel0 %>% length
    tnoise[sel0] -> t[sel0]
    t %>% matrix(nrow = nrow(Dat3),byrow=F) %>% 
      data.frame(rowID=Dat2 %>% dplyr::select(1) %>% unlist,.) %>% setNames(dat %>% colnames) -> Dato
  }
  # -----
  # without factor
  # -----
  if(is.null(value) & length(sel0)>0 & is.null(factor))
  {  
    # find non random 0 row
    Dat1 %>% t %>% data.frame %>% lapply("==",0) %>% lapply(which) %>% lapply(length) %>% unlist -> v0
    # minimum 0 by row to check randam disctribution
    if(is.null(ns)){ Dat0 %>% ncol %>% "%/%"(2) -> ns }
    v0 %>% "<="(ns) %>% which -> selkeep
    # chek randam 0 distribution across factor level
    c(selkeep) %>% unique -> selkeep2
    Dat1 %>% dplyr::slice(selkeep2) %>% data.frame(rowID=rowID[selkeep2],.) -> Dat2
    Dat2 %>% dplyr::select(-1) -> Dat3
    #
    # create background noise matrix
    Dat3 %>% unlist -> t
    t %>% "=="(0) %>% which -> sel00
    if(logbgvalue){ t[-sel00] %>% log2 %>% quantile(probs=quantv,na.rm=T) %>% as.numeric -> q100 }else
      { t[-sel00] %>% quantile(probs=quantv,na.rm=T) %>% as.numeric -> q100 }
    Dat3 %>% edgeR::DGEList() %>% edgeR::calcNormFactors(method = c("TMM")) -> dg0R
    Dat3 %>% t %>% edgeR::DGEList() %>% edgeR::calcNormFactors(method = c("TMMwsp")) -> dg0C
    dg0C$samples$norm.factors %o% dg0R$samples$norm.factors %>% "*"(q100) %>% data.frame -> dtnoise
    #
    # replace 0 by background noise
    Dat3 %>% unlist -> t
    dtnoise %>% unlist -> tnoise
    t %>% unlist %>% "=="(0) %>% which -> sel0
    tnoise[sel0] -> t[sel0]
    t %>% matrix(nrow = nrow(Dat3),byrow=F) %>% data.frame(rowID=Dat2 %>% dplyr::select(1) %>% unlist,.) %>% setNames(dat %>% colnames) -> Dato
  }
  Dato %>% return()
}
