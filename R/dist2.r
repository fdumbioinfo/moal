#' @title pearson distance dissimilarity
#' @description Function to compute the distance dissimilarity matrix using pearson correlation
#' @param x matrix or data.frame with numeric value
#' @param method character indicating distance method
#' @details Function to compute the distance dissimilarity matrix using pearson correlation
#'  distance methods available are "pearson", "kendall", "spearman" like in cor function (package stats)
#'  the matrix is firt transposate
#' @return matrix of dissimilarity (as.dist)
#' @examples
#' # not run
#' # dist2(rnorm(matrix(10*100, nrow = 10 )))
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @noRd
dist2 <- function( x , method = "pearson" )
{
  as.dist( abs( 1 - cor( t(x) , method = method ) )/2 )
}












