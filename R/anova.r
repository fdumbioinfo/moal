#' @title analysis of variance
#' @param dat numeric
#' @param model character
#' @param logratio logical change fc (by default) in log2ratio
#' @return data.frame
#' @examples
#' # not run
#' # library(magrittr)
#' # data(sif1)
#' # data(mat1)
#' # table(sif1$F1 , sif1$F2)
#' # mat1 %>% anova( sif1 , model = "F1+F2+F1*F2" )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm aov rnorm TukeyHSD
#' @importFrom utils capture.output
#' @importFrom broom tidy
#' @noRd
anova <- function( dat, model, logratio = FALSE )
{
  dat %>% names -> Names
  names(dat)[1] <- "y"
  # anova
  dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% tidy -> r0
  # pval
  r0 %>% select(.data$term,.data$p.value ) %>% slice(-nrow(.)) -> r1
  # stats
  r0$sumsq %>% "/"(.,sum(.)) %>% "*"(.,100) %>% setNames(paste("Sumsq_",r0$term %>% gsub(":","x",.),sep="")) -> Sumsq
  r0$meansq %>% "/"(.,sum(.)) %>% "*"(.,100) %>% setNames(paste("Meansq_",r0$term %>% gsub(":","x",.),sep="")) -> Meansq
  r0$statistic %>% replace(is.na(.),1) %>% setNames(paste("Fratio_",r0$term %>% gsub(":","x",.),sep="")) -> Fratio
  # contrasts
  dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% TukeyHSD %>% tidy %>% select(.data$contrast,.data$adj.p.value) -> c0
  # fc
  dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% TukeyHSD %>% tidy %>% select(.data$estimate) %>% "^"(2 ,.) %>% unlist -> ratio
  # dat %>% lm( formula = paste( "y~" , model ), data = .) %>% aov %>% TukeyHSD %>% tidy %>% select( .data$estimate ) %>% unlist -> ratio
  replace(ratio, which( ratio < 1 ), -1/ratio[ which( ratio < 1 ) ]) -> fc
  if(logratio){ ratio %>% log2 -> fc }
  # return
  r1 %>% select(.data$term) %>% unlist %>% as.character %>% gsub("-","vs",.) %>% gsub(":","x",.) %>% paste("p_",.,sep="") -> rNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-","vs",.) %>% gsub(":","-",.) %>% paste("p_",.,sep="") -> cNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-","vs",.) %>% gsub(":","-",.) %>% paste("fc_",.,sep="") -> fcNames
  r1 %>% select(.data$p.value) %>% unlist %>% setNames(rNames)  -> r2
  c( c0 %>% select(.data$adj.p.value) %>% unlist,fc) %>% setNames(c(cNames,fcNames)) -> cfc0
  cfc0 %>% length -> nbcol
  rep( c(1:(nbcol/2) ) , rep( 2 , nbcol/2 ) ) -> sel
  replace( sel , seq( 2 , nbcol , 2 ) , c(1:(nbcol/2))+(nbcol/2) ) -> sel
  cfc0[sel] -> cfc1
  list( c( r2 , cfc1 ), c( Sumsq , Meansq , Fratio ) ) %>% return()
}