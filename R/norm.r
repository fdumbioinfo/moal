#' @title Normalization function
#' @description normalization and log2
#' @param dat data.frame
#' @param method character apply quantile normalization by default see details for all method
#' @param log logical apply log base 2
#' @return data.frame
#' @details
#' Normalization method available: "TMM","TMMwsp","RLE","upperquartile" from EdgeR and "quantile","scale","cyclicloess"
#' 
#' @examples
#' # not run
#' # norm(dt)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom limma normalizeBetweenArrays voom
#' @importFrom edgeR DGEList normLibSizes
#' @export
norm <- function( dat = NULL, method = "quantile", log = TRUE )
{
  # 
  if(method == "none"){ dat -> dat1 }else{ dat -> dat1 }
  # limma method
  if(method %>% "=="(c("quantile","scale","cyclicloess")) %>% any)
  {
    dat %>% limma::normalizeBetweenArrays(method=method) %>% as.data.frame -> dat1
  }
  # EdgeR method
  if(method %>% "=="(c("TMM","TMMwsp","RLE","upperquartile")) %>% any)
  {
    dat %>% edgeR::DGEList() %>% edgeR::normLibSizes(method = c("TMM")) -> dg0
    dat %>% sweep(2,dg0$samples$norm.factors,"*") -> dat1
  }
  # log
  # if(log){ dat1 %>% "+"(.,1.01) %>% log2 -> dato }else{ dat1 -> dato }
  if(log)
  { 
    dat1 %>% unlist -> t 
    t %>% "<="(1) %>% which -> sel
    if(length(sel)>0)
    { 
      t[sel] <- 1.01
      t %>% matrix(nrow=nrow(dat1),byrow=F) %>% data.frame -> dato
    }
  }else{ dat1 -> dato }
  #
  dato %>% return()
}