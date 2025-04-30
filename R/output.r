#' @title export data.frame in tab file
#' @param dt data.frame
#' @param filename character
#' @examples
#' # not run
#' # output( dt )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom utils write.table
#' @export
output <- function( dt , filename )
{ write.table( x = dt , file = filename , sep = "\t", row.names = F, col.names = T, quote = F ) }