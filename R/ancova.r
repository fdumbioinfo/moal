#' @title analysis of variance
#' @param dat numeric
#' @param model character
#' @details directory with anova plot and results
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
#' @importFrom dplyr select slice
#' @importFrom rlang .data
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm aov rnorm TukeyHSD
#' @importFrom utils capture.output
#' @importFrom broom tidy
#' @noRd
ancova <- function( dat, model, contrast )
{
  dat %>% names -> Names
  names(dat)[1] <- "y"
  # anova
  dat %>% lm(formula=paste( "y~", model), data=.) %>% aov %>% broom::tidy(.) %>% data.frame -> r0
  # pval
  r0 %>% dplyr::select(.data$term, .data$p.value) %>% dplyr::slice(-nrow(.)) -> r1
  # stats
  r0$sumsq %>% "/"(.,sum(.)) %>% "*"(.,100) %>% setNames(paste("Sumsq_",r0$term %>% gsub(":","x",.),sep="")) -> Sumsq
  r0$meansq %>% "/"(.,sum(.)) %>% "*"(.,100) %>% setNames(paste("Meansq_",r0$term %>% gsub(":","x",.),sep=""))-> Meansq
  r0$statistic %>% replace(is.na(.),1) %>% setNames(paste("Fratio_",r0$term %>% gsub(":","x",.),sep="")) -> Fratio
  # contrasts
  dat %>% lm(formula = paste("y~",model), data=.) %>% aov %>% TukeyHSD(which=contrast) %>%
    broom::tidy(.) %>% dplyr::select(.data$contrast, .data$adj.p.value) -> c0
  # fc
  dat %>% lm(formula=paste("y~",model), data=.) %>% aov %>% TukeyHSD(which=contrast) %>%
    broom::tidy(.) %>% dplyr::select(.data$estimate) %>% "^"(2,.) %>% unlist -> ratio
  replace(ratio, which(ratio < 1), -1/ratio[which(ratio < 1)]) -> fc
  # return
  r1 %>% dplyr::select(.data$term) %>% unlist %>% as.character %>%
    gsub("-","vs",.) %>% gsub(":","x",.) %>% paste("p_",.,sep="") -> rNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-" ,"vs",.) %>% gsub(":","-",.) %>% paste("p_",.,sep="") -> cNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-","vs",.) %>% gsub(":","-",.) %>% paste("fc_",.,sep="") -> fcNames
  r1 %>% select( .data$p.value) %>% unlist %>% setNames( rNames  ) -> r2
  c(c0 %>% select( .data$adj.p.value ) %>% unlist, fc) %>% setNames(c( cNames, fcNames)) -> cfc0
  cfc0 %>% length -> nbcol
  rep(c(1:(nbcol/2)), rep(2,nbcol/2)) -> sel
  replace(sel, seq(2,nbcol,2), c(1:(nbcol/2))+(nbcol/2)) -> sel
  cfc0[sel] -> cfc1
  list(c(r2, cfc1), c(Sumsq, Meansq, Fratio)) %>% return()
}
