#' @title Hierarchical clustering classification
#' @description make a hierarchical clustering classification
#' @param dat matrix numeric
#' @param factor factor
#' @param title character
#' @param plot logical
#' @param method character to choose agglomerative method of clustering dendrogramm
#' @param legendtitle character
#' @param cexlabel numeric
#' @param legend logical display factor legend
#' @details
#'  Possible agglomerative method are the same hclust fonction : "complete" method by default
#' @return a dendrogram
#' @examples
#' # not run
#' # hc(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dendextend rotate set
#' @importFrom stats hclust as.dendrogram
#' @importFrom graphics plot legend
#' @noRd
hc <- function( dat, factor = NULL, title = "Hierarchical Clustering" , plot = TRUE,
                method = "complete", legendtitle = "TREATMENT", cexlabel = 0.45, legend = TRUE)
{
  dat %>% t -> dat
  dat %>% dist2 %>% stats::hclust(method=method) %>% stats::as.dendrogram(.) -> dend0
  #
  if( !is.null(factor) & plot & dim(dat)[1] == length(factor) )
  {
    dend0 %>% labels %>% match(rownames(dat)) %>% "["(factor,.) %>% factortocolor -> factorcol
    dend0 %>% dendextend::set("labels_colors", value=factorcol) %>% 
      dendextend::set("labels_cex",cexlabel) %>% plot(main=title,cex.main=0.8,ylab="distance correlation")
    # legend
    if(legend)
    {
      graphics::legend(
        "topright",title=legendtitle,legend=levels(factor),
        col=palette0[1:length(levels(factor))],lty=1,lwd=5,cex=0.7)
    }
    FALSE -> plot
  }
  #
  if(plot){ dend0 %>% dendextend::set("labels_cex",cexlabel ) %>% plot(main=title,cex.main=0.8,ylab="distance correlation") }
  dend0 %>% return()
}