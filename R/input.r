#' @title  import tab file in data.frame
#' @param filename character path to the file to read
#' @param sep character for field separator
#' @param quote character for field quote
#' @details wrapper of read.table function for tabular separated files
#' @return data.frame
#' @examples
#' # not run
#' # input( "filename" ) -> dt
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @export
input <- function( filename , sep = "\t" , quote = "" )
{
  read.table(
    filename, header=TRUE, sep=sep, quote=quote,
    check.names=FALSE, comment.char="",
    stringsAsFactors=FALSE ) %>% return()
}