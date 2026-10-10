#' @title analysis of variance
#' @param dat numeric
#' @param model character
#' @param logratio logical change fc (by default) in log2ratio
#' @param anovafc character fold-change calculate with anova model (by default) use meanfc or meanfcdelog
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
#' @importFrom broom tidy
#' @importFrom stats aggregate
#' @noRd
anova <- function( dat = NULL, model = NULL, logratio = FALSE, anovafc = "anovafc")
{
  i=1
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
  # dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% TukeyHSD %>% tidy %>% select(.data$estimate) %>% "^"(2 ,.) %>% unlist -> ratio
  dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% TukeyHSD %>% tidy %>% select(.data$estimate) %>% "^"(2 ,.) %>% unlist -> ratio
  dat %>% lm(formula=paste("y~",model),data=.) %>% aov %>% TukeyHSD %>% tidy %>% select(.data$estimate) %>% unlist -> Est0
  # dat %>% lm( formula = paste( "y~" , model ), data = .) %>% aov %>% TukeyHSD %>% tidy %>% select( .data$estimate ) %>% unlist -> ratio
  # replace(ratio, which( ratio < 1 ), -1/ratio[ which( ratio < 1 ) ]) -> fc
  # ratio < 1 equivalent to Est0 < 0 because 2 power 0 = 1
  # replace(ratio,which(Est0 < 0),-1/ratio[which(ratio < 1)]) -> fc
  replace(ratio,which(Est0 < 0),-1/ratio[which(Est0 < 0)]) -> fc
  # if(logratio){ ratio %>% log2 -> fc }
  if(logratio){ Est0 %>% log2 -> fc }
  # return
  r1 %>% select(.data$term) %>% unlist %>% as.character %>% gsub("-","vs",.) %>% gsub(":","x",.) %>% paste("p_",.,sep="") -> rNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-","vs",.) %>% gsub(":","-",.) %>% paste("p_",.,sep="") -> cNames
  c0 %>% select(.data$contrast) %>% unlist %>% gsub("-","vs",.) %>% gsub(":","-",.) %>% paste("fc_",.,sep="") -> fcNames
  r1 %>% select(.data$p.value) %>% unlist %>% setNames(rNames)  -> r2
  c(c0 %>% select(.data$adj.p.value) %>% unlist,fc) %>% setNames(c(cNames,fcNames)) -> cfc0
  cfc0 %>% length -> nbcol
  rep(c(1:(nbcol/2)),rep(2,nbcol/2)) -> sel
  replace(sel,seq(2,nbcol,2),c(1:(nbcol/2))+(nbcol/2)) -> sel
  cfc0[sel] -> cfc1
  #
  # fc without using anova estimate
  #
  # meanfc
  if(anovafc == "meanfc")
  {
    model %>% strsplit("\\+") %>% unlist %>% strsplit("\\*") %>% unlist %>% unique -> f
    f %>% paste0("^",.,"$") %>% paste0(collapse = "|") -> GrepF
    dat %>% colnames %>% grep(GrepF,.) %>% select(dat,.) -> ff
    c0 %>% select(.data$contrast) %>% unlist -> Contr0
    meanfc <- foreach(i=1:length(Contr0)) %do%
      {
        # Numerator
        Contr0[i] %>% strsplit("-") %>% unlist %>% "["(1) %>% strsplit(":") %>% unlist %>% lapply("==",ff) %>% lapply(which,arr.ind=T) %>% lapply("[",,1) %>% unlist -> fff
        if(fff %>% table %>% ">"(1) %>% which %>% any){ fff %>% table %>% ">"(1) %>% which %>% dat$y[.] %>% mean -> MeanN }else{ fff %>% dat$y[.] %>% mean -> MeanN }
        # Denominator
        Contr0[i] %>% strsplit("-") %>% unlist %>% "["(2) %>% strsplit(":") %>% unlist %>% lapply("==",ff) %>% lapply(which,arr.ind=T) %>% lapply("[",,1) %>% unlist -> fff
        if(fff %>% table %>% ">"(1) %>% which %>% any){ fff %>% table %>% ">"(1) %>% which %>% dat$y[.] %>% mean -> MeanD }else{ fff %>% dat$y[.] %>% mean -> MeanD }
        # fc
        ifelse((MeanN/MeanD) %>% ">"(1),MeanN/MeanD,-1/(MeanN/MeanD))
      }
    cfc1 %>% names %>% grep("^fc_",.) -> selfc
    meanfc %>% unlist -> cfc1[selfc] 
  }
  #
  # meanfcdelog
  if(anovafc == "meanfcdelog")
  {
    model %>% strsplit("\\+") %>% unlist %>% strsplit("\\*") %>% unlist %>% unique -> f
    f %>% paste0("^",.,"$") %>% paste0(collapse = "|") -> GrepF
    dat %>% colnames %>% grep(GrepF,.) %>% select(dat,.) -> ff
    c0 %>% select(.data$contrast) %>% unlist -> Contr0
    meanfc <- foreach(i=1:length(Contr0)) %do%
      {
        # Numerator
        Contr0[i] %>% strsplit("-") %>% unlist %>% "["(1) %>% strsplit(":") %>% unlist %>% lapply("==",ff) %>% lapply(which,arr.ind=T) %>% lapply("[",,1) %>% unlist -> fff
        if(fff %>% table %>% ">"(1) %>% which %>% any){ fff %>% table %>% ">"(1) %>% which %>% dat$y[.] %>% "^"(2,.) %>% mean -> MeanN }else{ fff %>% dat$y[.] %>% "^"(2,.) %>% mean -> MeanN }
        # Denominator
        Contr0[i] %>% strsplit("-") %>% unlist %>% "["(2) %>% strsplit(":") %>% unlist %>% lapply("==",ff) %>% lapply(which,arr.ind=T) %>% lapply("[",,1) %>% unlist -> fff
        if(fff %>% table %>% ">"(1) %>% which %>% any){ fff %>% table %>% ">"(1) %>% which %>% dat$y[.] %>% "^"(2,.) %>% mean -> MeanD }else{ fff %>% dat$y[.] %>% "^"(2,.) %>% mean -> MeanD }
        # fc
        ifelse((MeanN/MeanD) %>% ">"(1),MeanN/MeanD,-1/(MeanN/MeanD))
      }
    cfc1 %>% names %>% grep("^fc_",.) -> selfc
    meanfc %>% unlist -> cfc1[selfc]
  }
  #
  list(c(r2,cfc1),c(Sumsq,Meansq,Fratio)) %>% return()
}
