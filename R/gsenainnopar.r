#' @title Gene set enrichment analysis and interaction network
#' @description Gene set enrichment analysis and interaction network
#' @param omicdata character data.frame see details
#' @param keywords data.frame keywords list for geneset selection
#' @param species character hs mm rn dr ss
#' @param filtergeneset character list to filter MSigDB geneset collection
#' @param threshold numeric pval 0.05 fc 1.5 by default see details
#' @param dat data.frame file paths
#' @param factor factor factor for heatmap color
#' @param topdeg numeric top DEGs number to plot on network
#' @param rangedeg numeric top DEGs from 1 to topdeg by rangedeg to plot on network
#' @param topgeneset numeric top geneset number to plot on network
#' @param topena numeric top geneset for ena plot
#' @param twotailena logical do
#' @param bg numeric background used for functional analysis over-representation test
#' @param intmaxdh numeric maximum number of interaction to use for Davidson and Harel algorithm layout
#' @param nodesize numeric change Symbol size
#' @param doena logical do MSigDB enrichment analysis
#' @param gsearank character to choose gsea rank type among fc (by default) logration logfc sqrt 
#' @param layout numeric for layout neetwork 1 fr by default 2 dh 3 tree 4 circle 5 grid 6 sphere
#' @param dotopnetwork logical do top networks
#' @param dotopheatmap logical do top heatmap
#' @param dotopgenesetnetwork logical do geneset networks
#' @param dotopgenesetheatmap logical do geneset heatmap
#' @param dogmtgenesetnetwork logical do keyword networks
#' @param dogmtgenesetheatmap logical do keyword heatmap
#' @param path character for relative path of output directory
#' @param dirname character name for output
#' @param dopar logical TRUE for parallelization
#' @return file with enrichment analysis results
#' @details
#' 
#' omicdata needs a data.frame with at list 4 column: rowID, (p-values,fold-change) x N and Symbol annotation.
#' 
#' Symbol list are accepted to make ORA enrichment analysis.
#' 
#' To generate heatmap dat and factor parameter are needed. dat accepted complete matrix with rowID for first column. 
#' 
#' dat row IDs must match with omicdata row IDs. 
#' 
#' Make MSigDB enrichment analysis using GSEA method for non filtering list as input (> 2000)
#' 
#' Make MSigDB Over-Representation enrichment analysis (ORA) using Fisher exact test for list < 2000
#' 
#' Generate STRINGDB interaction network and heatmap for top geneset according to topena par (80 by default)
#' 
#' Only features with p-values < 0.05 et fold-change > 1.1 are displayed on geneset heatmaps (threshold = 1 by default).
#' 
#' See omic function details to display all threshold
#' 
#' @examples
#' # not run
#' # gsenainnopar( omicdata , species = "mm")
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange select desc last_col
#' @importFrom rlang .data
#' @importFrom stats fisher.test setNames p.adjust
#' @importFrom foreach foreach %do% %:%
#' @importFrom utils capture.output
#' @importFrom fgsea fgsea
#' @importFrom forcats fct_reorder
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @noRd
gsenainnopar <- function(
    omicdata = NULL, keywords = NULL, species = "hs", dat = NULL, factor = NULL,
    filtergeneset = NULL, threshold = 1 , topdeg = 80, rangedeg = NULL, topena = 80, twotailena = TRUE, 
    topgeneset = 80, intmaxdh = 5000, nodesize = 0.60, bg = 25000,
    doena = TRUE, gsearank = "logfc", layout = 1,
    dotopnetwork = TRUE, dotopgenesetnetwork = FALSE, dogmtgenesetnetwork = TRUE,
    dotopheatmap = TRUE, dotopgenesetheatmap = FALSE, dogmtgenesetheatmap = TRUE,
    path = NULL, dirname = NULL, dopar = TRUE)
{
  i=j=k=l=1
  ifelse(is.null(dirname),"ena" -> DirName0,paste("ena_",dirname,sep="") -> DirName0)
  if(is.null(path)){"." -> path}
  path %>% file.path(DirName0) -> Path0
  Path0 %>% dir.create
  file.path(Path0,"data") %>% dir.create
  if(is.null(dat)){ dotopheatmap <- F ; dotopgenesetheatmap <- F ; dogmtgenesetheatmap <- F}
  if(!is.null(threshold))
  {
    THRESHOLD <- TRUE
    thresholdlist[threshold] -> Threshold0
  }else{ THRESHOLD <- FALSE}
  # ----
  # 1- Annotation
  # ----
  ### example and avoid error with no data
  if(is.null(omicdata))
  {
    moalannotgene::genedbhs -> t
    t$Symbol %>% grep("^LOC",.,invert = T) %>% t$Symbol[.] -> tt
    tt %>% unique %>% sample(500) -> Symbol
    paste("row",1:500,sep="") -> rowID
    rep(0.05,500) -> p_AvsB
    rep(1.05,500) -> fc_AvsB
    data.frame(rowID,p_AvsB,fc_AvsB,Symbol=Symbol) -> omicdata
  }
  ### check for analysis without p-values or fold-changes.
  if(is.character(omicdata))
  {
    omicdata %>% as.character %>% unique -> Symbol
    paste("row",1:length(Symbol),sep="") -> rowID
    rep(0.05,length(Symbol)) -> p_AvsB
    rep(1.05,length(Symbol)) -> fc_AvsB
    data.frame(rowID,p_AvsB,fc_AvsB,Symbol=Symbol) -> omicdata
  }
  #
  mode(omicdata[,1]) <- "character"
  ### check Symbol
  omicdata$Symbol %>% as.character -> omicdata$Symbol 
  omicdata$Symbol %>% is.na %>% which -> sel
  if(length(sel) > 0){ omicdata %>% dplyr::slice(-sel) -> omicdata }
  omicdata$Symbol %>% grep("^$",.) -> sel
  if(length(sel) > 0){ omicdata %>% dplyr::slice(-sel) -> omicdata }
  omicdata$Symbol %>% "=="("NA") %>% which -> sel
  if(length(sel) > 0){ omicdata %>% dplyr::slice(-sel) -> omicdata }
  ### check NA p-value
  omicdata %>% colnames %>% grep("^p_",.) -> sel
  if(length(sel)>0)
  {
    foreach(i=1:length(sel)) %do%
      {
        omicdata[,sel] %>% unlist -> t
        t %>% is.na %>% which -> selNA
        if(length(selNA)){ t %>% replace(selNA,1) -> omicdata[,sel] }
      }
  }
  ### check NA fold-change
  omicdata %>% colnames %>% grep("^fc_",.) -> sel
  if(length(sel)>0)
  {
    foreach(i=1:length(sel)) %do%
      {
        omicdata[,sel] %>% unlist -> t
        t %>% is.na %>% which -> selNA
        if(length(selNA)){ t %>% replace(selNA,1) -> omicdata[,sel] }
      }
  }
  ### Annotation
  omicdata -> OmicData0
  OmicData0$Symbol %>% annot(dboutput="ebi",species=species,idtype = "SYMBOL") -> a0
  a0$Symbol -> OmicData0$Symbol
  # ebi annotation to add same ebi symbol exact match for network
  OmicData0 %>% data.frame(GeneID = a0$NCBIGeneID, check.names=F) -> OmicData1
  # ncbi geneid annotation to not loose GeneID not present in ebi
  OmicData0$Symbol %>% annot(dboutput="ncbi",species=species,idtype = "SYMBOL") -> a0GeneID
  OmicData1$GeneID <- a0GeneID$GeneID
  OmicData1$GeneID %>% is.na %>% "!"(.) %>% which -> sel
  if(length(sel) > 0){ OmicData1 %>% dplyr::slice(sel) -> OmicData2 }else{ OmicData1 -> OmicData2 } 
  OmicData2 %>% colnames %>% grep("p_",.) -> sel1
  t5 <- foreach(i=1:length(sel1)) %do%
    {
      OmicData2 %>% dplyr::select(1,c(sel1[i],sel1[i]+1),.data$Symbol,.data$GeneID ) -> t
      t %>% dplyr::arrange(.data[[colnames(t)[2]]]) -> tt
      tt %>% dplyr::group_by(.data$GeneID) %>% dplyr::slice(1) %>% data.frame(check.names=F) -> ttt
      ttt %>% dplyr::arrange(.data[[colnames(ttt)[2]]]) -> t4
      t4 %>% dplyr::group_by(.data$Symbol) %>% dplyr::slice(1) %>% data.frame(check.names=F)
    }
  # output
  t5 %>% lapply("[",1) %>% unlist %>% unique %>% data.frame(rowID=.) %>% inner_join(OmicData2) -> OutPut0
  OutPut0$Symbol %>% annot(dboutput = "ebi",species=species) -> a0
  a0 %>% dplyr::select(-.data$NCBIGeneID,-.data$Symbol) %>% data.frame(OutPut0,.,check.names=F) -> OutPut1
  paste("omicdata","_",species,"_",nrow(OutPut1),".tsv",sep="") -> FileName0
  OutPut1 %>% output(file.path(Path0,"data",FileName0))
  # foldchange
  foldchange <- foreach(i=1:length(t5)) %do%
    { t5[[i]] %>% dplyr::arrange(.data[[colnames(t5[[i]])[3]]]) %>% dplyr::select(4,3) }
  # pval
  pval <- foreach(i=1:length(t5)) %do%
    { t5[[i]] %>% dplyr::arrange(.data[[colnames(t5[[i]])[2]]]) %>% dplyr::select(4,2) }
  # ----
  # 2 - Top networks
  # ----
  if(dotopnetwork)
  {
    "NetworksDEGs" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    #
    foreach(j=1:length(t5)) %do%
      {
        t5[[j]] %>% dplyr::arrange(.data[[colnames(t5[[j]])[2]]]) %>% dplyr::slice(1:topdeg) -> p
        # output
        p$rowID %>% data.frame(rowID=.) %>% inner_join(OutPut1) -> Deg1
        paste("DEGs_top",topdeg,"_",colnames(t5[[j]])[2] %>% gsub("p_","",.),"_",nrow(Deg1),".tsv",sep="") -> FileName0
        # network
        Deg1 %>% output(file.path(Path0,DirName1,FileName0))
        if(nrow(Deg1) > 0)
        {
          Deg1$Symbol %>% sort %>% unique -> nodelist
          file.path(Path0,DirName1) -> pathn
          paste("top",topdeg,colnames(t5[[j]])[2] %>% gsub("p_","",.),sep="") -> filename
          networkgena4pval4(nodelist=nodelist,foldchange=foldchange,pval=pval,species=species,
                            intmaxdh=intmaxdh,filename=filename,path=pathn, layout = layout,
                            nodelabelsize=0.7,nodesize=0.7,edgeweight=0.1,edgewidth=0.7)
        }
      }
    #
    file.path(Path0,DirName1) %>% list.files(full.names = T) -> lf0
    if(length(lf0)>0)
    {
      lf0 %>% basename %>% grep("^DEGs_.*.tsv",.) %>% lf0[.] -> lf1
      file.copy(lf1,file.path(file.path(Path0,"data")))
      lf1 %>% file.remove
    }
    #
    file.path(Path0,DirName1) %>% list.files(full.names = T) -> lf0
    if(length(lf0)>0)
    {
      lf0 %>% basename %>% grep("^Edges_.*.tsv",.) %>% lf0[.] -> lf1
      IntAll0 <- foreach(l=1:length(lf1),.combine = "rbind") %do% { lf1[l] %>% input(.) }
      IntAll0 %>% unique -> IntAll1
      paste("Interactions_",DirName1,"_",nrow(IntAll1),".tsv",sep="") -> FileNameInt0
      IntAll1 %>% output(file.path(file.path(Path0,"data",FileNameInt0)))
      lf1 %>% file.remove
    } 
  }
  # ----
  # 3 - top heatmap
  # ----
  if(dotopheatmap)
  {
    "HeatmapDEGs" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    #
    foreach(j=1:length(t5)) %do%
      {
        t5[[j]] %>% dplyr::arrange(.data[[colnames(t5[[j]])[2]]]) %>% dplyr::slice(1:topdeg) -> p
        p %>% dplyr::select(1,4) %>% inner_join(dat) -> t
        t %>% head
        t5[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\1",.) -> GrN
        t5[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\2",.) -> GrD
        c(GrD,GrN) %>% paste("^",.,"$",sep="") %>% paste0(collapse = "|") -> grep
        factor %>% grep(grep,.) -> sel
        t %>% dplyr::select(c(sel+2)) -> tt
        tt %>% head
        factor[sel] -> Factor0
        paste("Heatmap_top",nrow(tt),"_",colnames(t5[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),sep="") -> Title0
        paste("Heatmap_top",nrow(tt),"_",colnames(t5[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
        file.path(Path0,DirName1,FileName0) -> FileName1
        # Symbol size
        5000 -> maxheatmap
        c(5,10,20,30,40,50,100,200,1000,maxheatmap) -> NRowheatmap0
        c(1,0.9,0.8,0.7,0.6,0.5,0.35,0.25,0.2,0.1) -> CexRow0
        if(nrow(tt) < maxheatmap)
        {nrow(tt) %>% "<="(NRowheatmap0) %>% which %>% min %>% CexRow0[.] -> CexRow1}else{ CexRow0[10] -> CexRow1 }
        #
        pdf(FileName1)
        tt %>% heatmap(Factor0,labCol=colnames(tt),labRow=t$Symbol,
                       dendrogram="row",cexRow=CexRow1,cexCol=0.5)
        title(main=Title0,cex.main=0.7)
        graphics.off()
      }
  }
  # ----
  # 4 - Ortholog annotation for MSigDB enrichment analysis and top geneset network
  # ----
  if(!(species == "hs"))
  {
    omicdata -> OmicData0
    OmicData0$Symbol %>% annot(species=species,dboutput="ncbi") -> a0
    a0$Symbol -> OmicData0$Symbol
    OmicData0 %>% data.frame(GeneID = a0$GeneID, check.names=F) -> OmicData1
    OmicData1$GeneID %>% is.na %>% which -> sel
    if(length(sel) > 0){ OmicData1 %>% dplyr::slice(-sel) -> OmicData2 }else{ OmicData1 -> OmicData2 } 
    OmicData2 %>% colnames %>% grep("p_",.) -> sel1
    t5 <- foreach(i=1:length(sel1)) %do%
      {
        OmicData2 %>% dplyr::select(1,c(sel1[i],sel1[i]+1),.data$Symbol,.data$GeneID ) -> t
        t %>% dplyr::arrange(.data[[colnames(t)[2]]]) -> tt
        tt %>% dplyr::group_by(.data$GeneID) %>% dplyr::slice(1) %>% data.frame(check.names=F) -> ttt
        ttt %>% dplyr::arrange(.data[[colnames(ttt)[2]]]) -> t4
        t4 %>% dplyr::group_by(.data$Symbol) %>% dplyr::slice(1) %>% data.frame(check.names=F)
      }
    # ortho
    t5[[1]] %>% dplyr::select(.data$Symbol) %>% unlist %>% annot(species=species,ortholog=T) -> a0
    a0 %>% dplyr::select(.data$GeneID,.data$OtherGeneID,.data$Symbol) %>% setNames(c("HsGeneID","GeneID","Symbol")) -> a1
    t6 <- foreach(i=1:length(t5)) %do%
      {
        t5[[i]] %>% dplyr::select(-.data$Symbol) -> t
        a1 %>% dplyr::inner_join(t,.) -> tt
        colnames(tt)[c(5,4)] <- c("GeneID",paste(species,"GeneID",sep=""))
        tt
      }
    # output
    OutPut0 <- foreach(i=1:length(t6),.combine = "rbind") %do% { t6[[i]][,c(1,6,5,4)] }
    OutPut0 %>% unique -> OutPut1
    OutPut1 %>% inner_join(OmicData0 %>% dplyr::select(-.data$Symbol),.) -> OutPut2
    OutPut2$Symbol %>% annot(dboutput = "ebi") -> a0
    a0 %>% dplyr::select(-.data$NCBIGeneID,-.data$Symbol) %>% data.frame(OutPut2,.,check.names=F) -> OutPut3
    paste("omicdata_ortholog_hs_",nrow(OutPut3),".tsv",sep="") -> FileName0
    OutPut3 %>% output(file.path(Path0,"data",FileName0))
    # foldchange
    foldchange <- foreach(i=1:length(t6)) %do% 
      { t6[[i]] %>% dplyr::arrange(.data[[colnames(t5[[i]])[3]]]) %>% dplyr::select(6,3) }
    # pval
    pval <- foreach(i=1:length(t6)) %do% 
      { t6[[i]] %>% dplyr::arrange(.data[[colnames(t5[[i]])[2]]]) %>% dplyr::select(6,2) }
  }else{t5 -> t6}
  # ----
  # 5 - MSigDB enrichment analysis 
  # ----
  ### gsea enrichment
  if(dotopgenesetnetwork | dotopgenesetheatmap){ doena <- T }
  if(doena & nrow(OutPut1) > 2000)
  {
    file.path(Path0,"ena") %>% dir.create
    moalannotgene::genesetdb -> t
    # t %>% names
    #
    # preprocessing
    gseastats0 <- foreach(j=1:length(t6)) %do%
      {
        t6[[j]] -> Omicdataf0
        # pval rank
        Omicdataf0 %>% dplyr::arrange(.data[[colnames(Omicdataf0)[2]]]) -> Omicdataf1
        # replace 0 value 
        Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf0)[2]]]) %>% "=="(0) %>% which -> sel
        if(length(sel) > 0)
        {
          Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf0)[2]]]) %>% "=="(0) %>% which -> sel
          colnames(Omicdataf0)[2] %>% grep(colnames(Omicdataf0)) -> selCol
          Omicdataf1[sel,selCol] <- 25
        }
        # twotailena = T
        Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[3]]]) %>% unlist %>% sign -> SignFC0
        1 -> ValueFC0
        Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[2]]]) %>% unlist %>% log10 %>% "*"(-1) %>% "*"(SignFC0) %>% "*"(ValueFC0) -> Signlog10pval0
        # decrease fold-change outliers ponderation
        Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[2]]]) %>% unlist %>% "<"(0.9) %>% which -> sel01
        if(length(sel01) > 0)
        {
          if(gsearank == "fc"){ Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[3]]]) %>% unlist %>% abs -> ValueFC0 }
          if(gsearank == "logratio"){ Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[3]]]) %>% unlist %>% fctolog2ratio %>% abs -> ValueFC0 }
          if(gsearank == "logfc"){ Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[3]]]) %>% unlist %>% abs %>% "+"(1) %>% log2 -> ValueFC0 }
          if(gsearank == "sqrt"){ Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf1)[3]]]) %>% unlist %>% abs %>% sqrt -> ValueFC0 }
          # 
          Signlog10pval0[sel01] %>% "*"(ValueFC0[sel01]) -> Signlog10pval0[sel01]
        }
        Omicdataf1 %>% data.frame(Signlog10pval0) %>% dplyr::arrange(Signlog10pval0) -> Omicdataf2
        Omicdataf2$Signlog10pval0 -> stats
        names(stats) <- Omicdataf2$GeneID
        # c(stats,CompName0)
        list(colnames(Omicdataf0)[2] %>% sub("p_","",.),stats)
      }
    # fsgsea
    gseastats1 <- foreach(i=1:length(t)) %do%
      {
        t[[i]] %>% lapply("[",2) 
        t[[i]] %>% lapply("[",1) %>% unlist -> GeneSetName0
        t[[i]] %>% lapply("[",4) %>% unlist(recursive = F) -> GeneSetNameList0
        names(GeneSetNameList0) <- GeneSetName0
        foreach(j=1:length(gseastats0)) %do%
          {
            GeneSetNameList0
            gseastats0[[j]][2] %>% unlist
            list(gseastats0[[j]][1] %>% as.character,gseastats0[[j]][2] %>% unlist,
                 GeneSetNameList0,t %>% names %>% "["(i))
          }
      }
    gseastats1 %>% unlist(recursive = F) -> gseastats2
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(gseastats2),.packages=c("magrittr","dplyr","moal","foreach","fgsea","stringr","ggplot2")) %do%
      {
        set.seed(123679)
        fgsea0 <- fgsea::fgsea(pathways=gseastats2[[i]][[3]],stats=gseastats2[[i]][[2]],minSize=15,maxSize=500,scoreType="std",nproc=1)
        if(nrow(fgsea0)>0)
        {
          fgsea0 -> fgsea1
          fgsea1$leadingEdge %>% lapply(paste0,collapse="|") %>% unlist(recursive = F) -> fgsea1$leadingEdge
          fgsea0$leadingEdge %>% lapply(length) %>% unlist -> OverlapSize
          fgsea1 %>% dplyr::select(1,2,6,7,8) %>% data.frame(OverlapSize=OverlapSize) %>%
            dplyr::select(c(1:(last_col()-2),last_col(),last_col()-1)) -> fgsea2
          c("Name","pval","NES","GeneSetSize","OverlapSize","OverlapGeneIDList") -> SetNames0
          fgsea2 %>% dplyr::arrange(.data$pval) %>% data.frame %>% setNames(SetNames0) -> fgsea3
          fgsea3 -> fgseapval
        }else{fgseapval <- NULL}
        #
        if(!is.null(fgseapval))
        {
          spf <- foreach(k=1:nrow(fgseapval)) %do%
            {
              fgseapval$OverlapGeneIDList[k] %>% strsplit("\\|") %>% unlist %>% "match"(Omicdataf0$GeneID) %>% 
                Omicdataf0$Symbol[.] %>% paste0(collapse = "|")
            }
          fgseapval %>% data.frame(OverlapSymbolList=spf %>% unlist(recursive = F)) %>%
            dplyr::select(c(1,5,4,7,2,3,6)) -> fgseapval1
          gseastats2[[i]][[4]] %>% strsplit("\\|") %>% unlist %>% "["(2) %>% stringr::str_to_upper(.) -> DirName1
          paste(DirName1,"_ena_",gseastats2[[i]][[1]],"_",nrow(fgseapval1),".tsv",sep="") -> FileName0
          fgseapval1 %>% output(file.path(Path0,"ena",FileName0))
          # plot
          fgseapval1 %>% dplyr::slice(1:topena) -> fgseapval1plot
          if(fgseapval1plot %>% nrow %>% ">"(0))
          {
            fgseapval1plot$pval %>% log10 %>% "*"(-1) -> Log10Pval
            fgseapval1plot %>% data.frame(Log10Pval) -> fgseapval1plot2
            # reduce long geneset name
            fgseapval1plot2$Name %>% as.character %>% nchar -> Nchar0
            Nchar0 %>% ">"(50) %>% which -> selNchar
            if(length(selNchar)>0)
            {
              fgseapval1plot2$Name[selNchar] %>% substr(1,45) -> Head0
              fgseapval1plot2$Name[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
              fgseapval1plot2$Name %>% as.character -> NameNchar0
              fgseapval1plot2$Name[selNchar] <- paste(Head0,Tail0,sep="...")
            }
            fgseapval1plot2 %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,Log10Pval)) -> fgseapval1plot3
            # barplot 
            fgseapval1plot3 %>% ggplot( aes(x=Log10Pval,y=.data$Name,fill=.data$NES)) -> p
            p + geom_bar(stat="identity",width = 0.95) -> p
            p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
            # p + theme_minimal() -> p
            # p + theme_light() -> p
            # p + theme_linedraw() -> p
            p + theme_bw() -> p
            p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 5),
                      axis.text.y = element_text(angle = 0, vjust = 0.4, hjust=1,size = 6)) -> p
            p + labs(x="log10(p-value)",y="Genesets") -> p
            p + ggtitle(DirName1,subtitle=gseastats2[[i]][[1]]) -> p
            p + theme(axis.title.x=element_text(size=12,face="bold"),
                      axis.title.y=element_text(size=12),axis.text.y=element_text(size=6),
                      plot.title=element_text(size=8,hjust=0.5),plot.subtitle=element_text(size=8,hjust=0.5)) -> p
            p + geom_vline(xintercept=-log10(0.05),color="darkgoldenrod1",size=0.2,linetype="dashed") -> p
            # output plot
            paste(DirName1,"_ena_",gseastats2[[i]][[1]],"_",nrow(fgseapval1plot3),".pdf",sep="") -> FileName0
            ggsave(plot=p,filename=file.path(Path0,"ena",FileName0))
          }
        }
      }
    if(dopar){parallel::stopCluster(cl);doParallel::stopImplicitCluster()}
    # 
    Path0 %>% file.path("ena") %>% list.files(full.names = T) -> l0
    l0
    t %>% names
    if(length(l0) > 0)
    {
      Path0 %>% file.path("ena") %>% list.files(full.names = T) -> l0
      # pathways
      l0 %>% basename %>% grep("^REACTOME|^KEGG|^WIKI|^HALLMAR|^KEGG|^CGP|^PID|^BIOCARTA",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","PATHWAYS") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","PATHWAYS")) 
      }
      # Gene Ontology
      l0 %>% basename %>% grep("^GOBP|^GOCC|^GOMM|^HPO",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","GO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","GO")) 
      }
      # LOCALISATION
      l0 %>% basename %>% grep("^CHR|^CELLTYPE",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","LOCALIZATION") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","LOCALIZATION")) 
      }
      # ONCOPATTERN
      l0 %>% basename %>% grep("^ONCO|^CGN|^CM",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","ONCO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","ONCO")) 
      }
      # IMMUNOPATTERN
      l0 %>% basename %>% grep("^VAX|^IMMUNO",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","IMMUNO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","IMMUNO")) 
      }
      # REGULATORS
      l0 %>% basename %>% grep("^TFT|^MIR",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","REGULATORS") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","REGULATORS")) 
      }
    l0 %>% file.remove
    }
  }
  #
  # ORA enrichment
  #
  if(doena & nrow(OutPut1) < 2000)
  {
    file.path(Path0,"ena") %>% dir.create
    moalannotgene::genesetdb -> Genesetdb0
    Genesetdb0 %>% names
    #
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(Genesetdb0),.packages=c("magrittr","dplyr","moal","foreach","stringr","ggplot2")) %do%
      {
        # gene set
        Genesetdb0[[i]] %>% lapply("[",1) %>% unlist -> GeneSetName0
        Genesetdb0[[i]] %>% lapply("[",4) %>% unlist(recursive = F) -> GeneSetNameList0
        names(GeneSetNameList0) <- GeneSetName0
        # 
        foreach(j=1:length(t6)) %do%
          {
            t6[[j]] -> Omicdataf0
            # pval rank
            Omicdataf0 %>% dplyr::arrange(.data[[colnames(Omicdataf0)[2]]]) -> Omicdataf1
            # 0 value 
            Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf0)[2]]]) %>% "=="(0) %>% which -> sel
            if(length(sel) > 0)
            {
              Omicdataf1 %>% dplyr::select(.data[[colnames(Omicdataf0)[2]]]) %>% "=="(0) %>% which -> sel
              colnames(Omicdataf0)[2] %>% grep(colnames(Omicdataf0)) -> selCol
              Omicdataf1[sel,selCol] <- 25
            }
            #
            1 -> overlapmin
            1.1 -> enascoremin
            Omicdataf1$GeneID -> GeneidList0
            Omicdataf1$Symbol -> SymbolList0
            Ena0 <- foreach(k=1:length(Genesetdb0[[i]]),.combine="rbind") %do%
              {
                Genesetdb0[[i]][[k]][[1]] -> GenesetName
                GeneidList0 %in% Genesetdb0[[i]][[k]][[4]] %>% which %>% GeneidList0[.] -> OverlapGeneIDList0
                GeneidList0 %in% Genesetdb0[[i]][[k]][[4]] %>% which %>% SymbolList0[.] -> OverlapSymbolList0
                Omicdataf1$Symbol %>% length -> ListSize0
                OverlapSymbolList0 %>% length -> OverlapSize0
                Genesetdb0[[i]][[k]][[4]] %>% length -> GenesetSize0
                ( OverlapSize0 / GenesetSize0 ) %>% round(2) -> OverlapRatio0
                ( ( OverlapSize0 / ListSize0 ) / (  GenesetSize0 / bg ) ) %>% round(2) -> EnaScore0
                if(  OverlapSize0 >= overlapmin & EnaScore0 > enascoremin )
                {
                  c(OverlapSize0,ListSize0,GenesetSize0,bg) %>% matrix(ncol=2) -> ContTable0
                  ContTable0 %>% fisher.test %>% unlist %>% "["(1) %>% as.numeric -> Pval0
                  c(
                    GenesetName,
                    paste0( OverlapSymbolList0, collapse = "|" ),
                    paste0( OverlapGeneIDList0, collapse = "|" ),
                    OverlapSize0,GenesetSize0,OverlapRatio0,EnaScore0,Pval0)
                }
              }
            #
            if(!is.null(Ena0))
            {
              Ena0 %>% data.frame(stringsAsFactors=F) %>%
              stats::setNames(c("Name","SymbolList","GeneIDList","OverlapSize","GenesetSize","OverlapRatio","ENAScore","pval")) -> Ena1
              Ena1 %>% colnames
              Ena1$OverlapSize %>% as.numeric -> Ena1$OverlapSize
              Ena1$GenesetSize %>% as.numeric -> Ena1$GenesetSize
              Ena1$OverlapRatio %>% as.numeric -> Ena1$OverlapRatio
              Ena1$ENAScore %>% as.numeric -> Ena1$ENAScore
              Ena1$pval %>% as.numeric -> Ena1$pval
              Ena1 %>% mutate( pvalFDR = .data$pval %>% p.adjust(method = "fdr") ) %>%
                mutate( log10pvalFDR = .data$pvalFDR %>% log10 %>% '*'( . , -1 ) %>% round( . , 4) ) %>%
                arrange( .data$pvalFDR  ) -> Ena1
              Ena1 %>% colnames  
              Ena1 %>% dplyr::select(c(1,4,5,2,8,7,3)) %>% 
                setNames(c("Name","OverlapSize","GeneSetSize","OverlapSymbolList","pval","NES","OverlapGeneIDList")) -> Ena1
              Ena1 %>% dplyr::arrange(.data$pval) -> Ena1
              Ena1 -> fgseapval1
              Genesetdb0 %>% names %>% "["(i) %>% strsplit("\\|") %>% unlist %>% "["(2) %>% stringr::str_to_upper(.) -> DirName1
              paste(DirName1,"_ena_",colnames(Omicdataf1)[3] %>% gsub("fc_","",.),"_",nrow(fgseapval1),".tsv",sep="") -> FileName0
              fgseapval1 %>% output(file.path(Path0,"ena",FileName0))
              #
              fgseapval1 %>% dplyr::slice(1:topena) -> fgseapval1plot
              if(fgseapval1plot %>% nrow %>% ">"(0))
              {
                fgseapval1plot$pval %>% log10 %>% "*"(-1) -> Log10Pval
                fgseapval1plot %>% data.frame(Log10Pval) -> fgseapval1plot2
                # reduce long geneset name
                fgseapval1plot2$Name %>% as.character %>% nchar -> Nchar0
                Nchar0 %>% ">"(50) %>% which -> selNchar
                if(length(selNchar)>0)
                {
                  fgseapval1plot2$Name[selNchar] %>% substr(1,45) -> Head0
                  fgseapval1plot2$Name[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
                  fgseapval1plot2$Name %>% as.character -> NameNchar0
                  NameNchar0
                  fgseapval1plot2$Name[selNchar] <- paste(Head0,Tail0,sep="...")
                }
                fgseapval1plot2 %>% dplyr::mutate(Name=forcats::fct_reorder(.data$Name,Log10Pval)) -> fgseapval1plot3
                # barplot 
                fgseapval1plot3 %>% ggplot( aes(x=Log10Pval,y=.data$Name,fill=.data$NES)) -> p
                p + geom_bar(stat="identity") -> p
                p + scale_color_gradient2(low="blue",mid="white",high="red",aesthetics="fill") -> p
                p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,size = 5),
                          axis.text.y = element_text(angle = 0, vjust = 1, hjust=1,size = 6)) -> p
                p + labs(x="log10(p-value)",y="Genesets") -> p
                p + theme(axis.title.x =element_text(size=12,face="bold"),axis.title.y=element_text(size=12))
                p + ggtitle(DirName1) -> p
                p + geom_vline(xintercept=-log10(0.05),color="orange",size=0.3) -> p
                p
                # output plot
                paste(DirName1,"_ena_",colnames(Omicdataf1)[3] %>% gsub("fc_","",.),"_",nrow(fgseapval1plot3),".pdf",sep="") -> FileName0
                ggsave(plot=p,filename=file.path(Path0,"ena",FileName0))
              }
            }
          }
      }
    if(dopar){parallel::stopCluster(cl);doParallel::stopImplicitCluster()}
    #
    Path0 %>% file.path("ena") %>% list.files(full.names = T,recursive = T) -> l0
    l0
    t %>% names
    if(length(l0) > 0)
    {
      # pathways
      l0 %>% basename %>% grep("^REACTOME|^KEGG|^WIKI|^HALLMAR|^KEGG|^CGP|^PID|^BIOCARTA",.) -> selgs
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","PATHWAYS") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","PATHWAYS")) 
      }
      # Gene Ontology
      l0 %>% basename %>% grep("^GOBP|^GOCC|^GOMM|^HPO",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","GO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","GO")) 
      }
      # LOCALISATION
      l0 %>% basename %>% grep("^CHR|^CELLTYPE",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","LOCALIZATION") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","LOCALIZATION")) 
      }
      # ONCOPATTERN
      l0 %>% basename %>% grep("^ONCO|^CGN|^CM",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","ONCO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","ONCO")) 
      }
      # IMMUNOPATTERN
      l0 %>% basename %>% grep("^VAX|^IMMUNO",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","IMMUNO") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","IMMUNO")) 
      }
      # REGULATORS
      l0 %>% basename %>% grep("^TFT|^MIR",.) -> selgs   
      if(length(selgs)>0)
      {
        Path0 %>% file.path("ena","REGULATORS") %>% dir.create
        l0[selgs] %>% file.copy(file.path(Path0,"ena","REGULATORS")) 
      }
      l0 %>% file.remove
    }
  }
  # ----
  # 6 - geneset networks from ena results
  # ----
  if(dotopgenesetnetwork)
  {
    moalannotgene::genesetdb -> Genesetdb0
    topgeneset -> Topgenesetval0
    if(!is.null(filtergeneset))
    {
      filtergeneset %>% paste0(collapse = "|") -> GrepFiltergeneset0
      Genesetdb0 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% unlist %>% 
        grep(GrepFiltergeneset0,.,ignore.case=T) %>% Genesetdb0[.] -> Genesetdb1
    }else{Genesetdb0 -> Genesetdb1}
    Genesetdb1 %>% names
    #
    "NetworkTopGenesets" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    GenesetNames <- foreach(k=1:length(Genesetdb1),.combine = "c") %do% { Genesetdb1[[k]] %>% lapply("[",1) %>% unlist }
    GenesetLists <- foreach(k=1:length(Genesetdb1),.combine = "c") %do% { Genesetdb1[[k]] %>% lapply("[",4) }
    Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% unlist %>%
      rep(Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",3) %>% unlist %>% as.numeric) -> GenesetCollection
    #
    file.path(Path0,"ena") %>% list.files(full.names = T,recursive = T) -> l0
    l0 %>% basename %>% grep("^.*_ena_.*.tsv$",.,value = F) %>% l0[.]  -> l1
    Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% stringr::str_to_upper(.) -> Collection0
    file.path(Path0,DirName1,Collection0) %>% lapply(dir.create)
    #
    TopGeneset0 <- foreach(i=1:length(l1),.combine = "c") %do%
      {
        l1[i] %>% input(.) -> ll0
        ll0 %>% dplyr::arrange()
        if(nrow(ll0) < topgeneset){ nrow(ll0) -> Topgenesetval0 }
        ll0$Name[1:Topgenesetval0] -> t
        topgeneset -> Topgenesetval0
        t
      }
    GenesetNames %>% "%in%"(TopGeneset0) %>% which -> sel
    GenesetLists[sel] %>% unlist(recursive = F) -> TopGeneset1
    GenesetNames[sel] %>% paste("|",GenesetCollection[sel] %>% stringr::str_to_upper(.) ,sep="") -> GenesetNames0
    names(TopGeneset1) <- GenesetNames0
    TopGeneset1 %>% lapply(length)
    ### threshold nodelist with pval < 0.05 & convert Geneset geneID in Symbol
    if(!THRESHOLD)
    {
      TopGeneset2 <- foreach(i=1:length(TopGeneset1),.combine = "c") %do%
        {
          TopGeneset1[i] -> t
          t %>% names
          t[[1]] %>% annot(.) -> tt
          t[[1]] <- tt$Symbol
          t
        }
    }
    if(THRESHOLD)
    {
      Pvalnodelist0 <- foreach(i=1:length(pval)) %do% 
        { 
          pval[[i]] %>% dplyr::arrange(.data[[colnames(pval[[i]])[2]]]) %>% 
            dplyr::filter(.data[[colnames(pval[[i]])[2]]] < unlist(Threshold0[[1]][1])) %>% dplyr::select(.data$Symbol) %>% unlist
        }
      Pvalnodelist0 %>% unlist %>% unique -> Pvalnodelist1
      TopGeneset2 <- foreach(i=1:length(TopGeneset1),.combine = "c") %do%
        {
          TopGeneset1[i] -> t
          t %>% names
          t[[1]] %>% annot(.) -> a0
          a0$Symbol %>% data.frame(Symbol=.) %>% inner_join(data.frame(Symbol=Pvalnodelist1)) -> tt
          tt$Symbol
          t[[1]] <- tt$Symbol
          t
        }
    }
    #
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(TopGeneset2),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %do%
      {
        TopGeneset2[[i]] %>% sort -> nodelist
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> titles
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameTop
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        file.path(Path0,DirName1,DirNameN0) -> pathn
        # intmaxdh <- 10000
        networkgena4pval4(nodelist=nodelist,foldchange=foldchange,pval=pval,species="hs",
                          intmaxdh=intmaxdh,filename=FileNameTop,path=pathn,title=title, layout=layout,
                          nodelabelsize=0.7,nodesize=0.7,edgeweight=0.1,edgewidth=0.7)
      }
    if(dopar){parallel::stopCluster(cl); doParallel::stopImplicitCluster()}
    #
    file.path(Path0,DirName1) %>% list.files(full.names = T, recursive = T) -> lf0
    if(length(lf0) > 0)
    {
      lf0 %>% basename %>% grep("^Edges_.*.tsv",.) %>% lf0[.] -> lf1
      IntAll0 <- foreach(l=1:length(lf1),.combine = "rbind") %do% { lf1[l] %>% input(.) }
      IntAll0 %>% unique -> IntAll1
      paste("Interactions_",DirName1,"_",nrow(IntAll1),".tsv",sep="") -> FileNameInt0
      IntAll1 %>% output(file.path(file.path(Path0,"data",FileNameInt0)))
      lf1 %>% file.remove
    }
  }
  # ----
  # 7 - geneset heatmap 
  # ----
  if(dotopgenesetheatmap)
  {
    moalannotgene::genesetdb -> Genesetdb0
    topgeneset -> Topgenesetval0
    if(!is.null(filtergeneset))
    {
      filtergeneset %>% paste0(collapse = "|") -> GrepFiltergeneset0
      Genesetdb0 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% unlist %>% 
        grep(GrepFiltergeneset0,.,ignore.case=T) %>% Genesetdb0[.] -> Genesetdb1
    }else{Genesetdb0 -> Genesetdb1}
    Genesetdb0 %>% names
    "HeatmapTopGenesets" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    GenesetNames <- foreach(k=1:length(Genesetdb1),.combine = "c") %do% { Genesetdb1[[k]] %>% lapply("[",1) %>% unlist }
    GenesetLists <- foreach(k=1:length(Genesetdb1),.combine = "c") %do% { Genesetdb1[[k]] %>% lapply("[",4) }
    Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% unlist %>%
      rep(Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",3) %>% unlist %>% as.numeric) -> GenesetCollection
    #
    file.path(Path0,"ena") %>% list.files(full.names=T,recursive=T) -> l0
    l0 %>% basename %>% grep("^.*_ena_.*.tsv$",.,value = F) %>% l0[.]  -> l1
    Genesetdb1 %>% names %>% strsplit("\\|") %>% lapply("[",2) %>% stringr::str_to_upper(.) -> Collection0
    file.path(Path0,DirName1,Collection0) %>% lapply(dir.create)
    TopGeneset0 <- foreach(i=1:length(l1),.combine = "c") %do%
      {
        l1[i] %>% input(.) -> ll0
        ll0 %>% dplyr::arrange()
        if(nrow(ll0) < topgeneset){ nrow(ll0) -> Topgenesetval0 }
        ll0$Name[1:Topgenesetval0] -> t
        topgeneset -> Topgenesetval0
        t
      }
    GenesetNames %>% "%in%"(TopGeneset0) %>% which -> sel
    GenesetLists[sel] %>% unlist(recursive = F) -> TopGeneset2
    GenesetNames[sel] %>% paste("|",GenesetCollection[sel] %>% stringr::str_to_upper(.) ,sep="") -> GenesetNames0
    names(TopGeneset2) <- GenesetNames0
    #
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(TopGeneset2),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %do%
      {
        TopGeneset2[[i]] %>% annot(.) -> a0
        a0$Symbol -> nodelist
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameTop
        # reduce long geneset name
        FileNameTop %>% as.character %>% nchar -> Nchar0
        Nchar0 %>% ">"(50) %>% which -> selNchar
        if(length(selNchar)>0)
        {
          FileNameTop[selNchar] %>% substr(1,45) -> Head0
          FileNameTop[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
          FileNameTop %>% as.character -> NameNchar0
          NameNchar0
          FileNameTop <- paste(Head0,Tail0,sep="...")
        }
        TopGeneset2[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        file.path(Path0,DirName1,DirNameN0) -> pathn
        foreach(j=1:length(t6)) %do%
          {
            # threshold
            if(THRESHOLD)
            {
              t6[[j]] %>% dplyr::arrange(.data[[colnames(t5[[j]])[2]]]) %>% 
                dplyr::filter((.data[[colnames(t5[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                                 .data[[colnames(t5[[j]])[3]]] > unlist(Threshold0[[1]][2]) )|(
                                   .data[[colnames(t5[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                                     .data[[colnames(t5[[j]])[3]]] < -unlist(Threshold0[[1]][2]))) -> p
              
            }else{ t6[[j]] -> p }
            #
            if(species=="hs"){p %>% dplyr::select(1,4) %>% inner_join(dat) -> t }else{p %>% dplyr::select(1,6) %>% inner_join(dat) -> t}
            nodelist %>% data.frame(Symbol=.) %>% inner_join(t) -> tt
            #
            if(nrow(tt) > 2)
            {
              t6[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\1",.) -> GrN
              t6[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\2",.) -> GrD
              c(GrD,GrN) %>% paste("^",.,"$",sep="") %>% paste0(collapse = "|") -> grep
              factor %>% grep(grep,.) -> sel
              tt %>% dplyr::select(c(sel+2)) -> ttt
              # factor[sel] %>% as.character %>% ordered(c(GrD,GrN)) -> Factor0
              factor[sel] -> Factor0
              #
              if(THRESHOLD)
              {
                paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",FileNameTop,"_",
                      colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
                paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",
                      colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),sep="") -> Title0
              }else
                {
                  paste("Heatmap_",FileNameTop,"_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
                  paste("Heatmap_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),sep="") -> Title0
                }
              Title0 %>% paste("\n\n",FileNameTop) -> Title1
              file.path(Path0,DirName1,DirNameN0,FileName0) -> FileName1
              # Symbol size
              5000 -> maxheatmap
              c(5,10,20,30,40,50,100,200,1000,maxheatmap) -> NRowheatmap0
              c(1,0.9,0.8,0.7,0.6,0.5,0.35,0.25,0.2,0.1) -> CexRow0
              if(nrow(ttt) < maxheatmap)
              {nrow(ttt) %>% "<="(NRowheatmap0) %>% which %>% min %>% CexRow0[.] -> CexRow1}else{ CexRow0[10] -> CexRow1 }
              #
              pdf(FileName1)
              ttt %>% heatmap(Factor0,labCol=colnames(ttt),labRow=tt$Symbol,
                              dendrogram="row",cexRow=CexRow1,cexCol=0.5)
              title(main=Title1,cex.main=0.6)
              graphics.off()
            }
          }
        }
    if(dopar){parallel::stopCluster(cl); doParallel::stopImplicitCluster()}
  }
  # ----
  # 8 - GMT network
  # ----
  if(!is.null(keywords) & dogmtgenesetnetwork)
  {
    "NetworksGMT" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    GenesetLists1 <- foreach(i=1:length(keywords), .combine = "c") %do%
      {
        keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
        file.path(Path0,DirName1,Keyword0) %>% dir.create
        keywords[i] %>% file("r") -> f 
        readLines(f,warn=FALSE) -> rl0
        close(f)
        rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2)) -> GenesetLists0
        rl0 %>% strsplit("\t") %>% lapply("[",1) %>% unlist -> Names0
        paste(Names0,"|",Keyword0,sep="") -> Names1
        names(GenesetLists0) <- Names1
      }
    #
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(GenesetLists1),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %do%
      {
        GenesetLists1[[i]] -> symbollist
        GenesetLists1[[i]] %>% annot(.,idtype = "SYMBOL") -> a0
        a0$Symbol -> nodelist
        GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameGeneset
        GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        file.path(Path0,DirName1,DirNameN0) -> pathn
        networkgena4pval4(nodelist=nodelist,foldchange=foldchange,pval=pval,species="hs",
                          intmaxdh=intmaxdh,filename=FileNameGeneset,path=pathn,title=title,layout = layout,
                          nodelabelsize=0.7,nodesize=0.7,edgeweight=0.1,edgewidth=0.7)
      }
    if(dopar){parallel::stopCluster(cl); doParallel::stopImplicitCluster()}
    #
    file.path(Path0,DirName1) %>% list.files(full.names=T,recursive=T) -> lf0
    lf0 %>% basename %>% grep("^Edges.*.tsv",.) %>% lf0[.] -> lf1
    IntAll0 <- foreach(l=1:length(lf1),.combine = "rbind") %do% { lf1[l] %>% input(.) }
    IntAll0 %>% unique -> IntAll1
    paste("Interactions_",DirName1,"_",nrow(IntAll1),".tsv",sep="") -> FileNameInt0
    IntAll1 %>% output(file.path(Path0,"data",FileNameInt0))
    lf1 %>% file.remove
  }
  # ----
  # 9 - GMT heatmap
  # ----
  if(!is.null(keywords) & dogmtgenesetheatmap)
  {
    "HeatmapGMT" -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    # moal:::thresholdlist[threshold] -> Threshold0
    keywords %>% basename %>% sub("(.*).gmt","\\1",.) -> DirNameGMT0
    file.path(Path0,DirName1,DirNameGMT0) %>% lapply(dir.create)
    # GenesetLists1 <- foreach(i=1:length(keywords), .combine = "c") %do%
    #   {
    #     keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
    #     keywords[i] %>% file("r") -> f
    #     readLines(f,warn=FALSE) -> rl0
    #     close(f)
    #     rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2)) -> GenesetLists0
    #     rl0 %>% strsplit("\t") %>% lapply("[",1) %>% unlist -> Names0
    #     paste(Names0,"|",Keyword0,sep="") -> Names1
    #     names(GenesetLists0) <- Names1
    #   }
    GenesetLists1 <- foreach(i=1:length(keywords), .combine = "c") %do%
      {
        keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
        keywords[i] %>% file("r") -> f
        readLines(f,warn=FALSE) -> rl0
        close(f)
        rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2)) -> GenesetLists0
        rl0 %>% strsplit("\t") %>% lapply("[",1) %>% unlist -> Names0
        paste(Names0,"|",Keyword0,sep="") -> Names1
        names(GenesetLists0) <- Names1
      }
    GenesetLists1 %>% length
    GenesetListsSymbol <- foreach(i=1:length(keywords), .combine = "c") %do%
      {
        keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
        keywords[i] %>% file("r") -> f
        readLines(f,warn=FALSE) -> rl0
        close(f)
        rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2))
      }
    GenesetListsSymbol %>% head
    GenesetListsSymbol %>% length
    GenesetListsSymbol %>% lapply(length)
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(GenesetLists1),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %do%
      {
        # GenesetLists1[[i]] -> symbollist
        GenesetListsSymbol[[i]] -> symbollist
        symbollist
        symbollist %>% length
        # GenesetLists1[[i]] %>% annot(.,idtype = "SYMBOL") -> a0
        symbollist %>% annot(.,idtype = "SYMBOL") -> a0
        a0$Symbol -> nodelist
        # GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameGeneset
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        # reduce long geneset name
        FileNameGeneset %>% as.character %>% nchar -> Nchar0
        Nchar0 %>% ">"(50) %>% which -> selNchar
        if(length(selNchar)>0)
        {
          FileNameGeneset[selNchar] %>% substr(1,45) -> Head0
          FileNameGeneset[selNchar] %>% substr(Nchar0-10,Nchar0)-> Tail0
          FileNameGeneset %>% as.character -> NameNchar0
          NameNchar0
          FileNameGeneset <- paste(Head0,Tail0,sep="...")
        }
        # GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        #
        foreach(j=1:length(t6)) %do%
          {
            # threshold
            if(THRESHOLD)
            {
              t6[[j]] %>% dplyr::arrange(.data[[colnames(t6[[j]])[2]]]) %>%
                dplyr::filter((.data[[colnames(t6[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                                 .data[[colnames(t6[[j]])[3]]] > unlist(Threshold0[[1]][2]) )|(
                                   .data[[colnames(t6[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                                     .data[[colnames(t5[[j]])[3]]] < -unlist(Threshold0[[1]][2]))) -> p
              
            }else{ t6[[j]] -> p }
            if(species=="hs"){p %>% dplyr::select(1,4) %>% inner_join(dat) -> t }else{p %>% dplyr::select(1,6) %>% inner_join(dat) -> t}
            nodelist %>% data.frame(Symbol=.) %>% inner_join(t) -> tt
            if(nrow(tt) > 2)
            {
              t6[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\1",.) -> GrN
              t6[[j]] %>% colnames %>% "["(2) %>% sub("^p_(.*)vs(.*)","\\2",.) -> GrD
              c(GrD,GrN) %>% paste("^",.,"$",sep="") %>% paste0(collapse = "|") -> grep
              factor %>% grep(grep,.) -> sel
              tt %>% dplyr::select(c(sel+2)) -> ttt
              # factor[sel] %>% as.character %>% ordered(c(GrD,GrN)) -> Factor0
              factor[sel] -> Factor0
              if(THRESHOLD)
              {
                paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",FileNameGeneset,"_",
                      colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
                # paste("Heatmap_",FileNameGeneset,"_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
                paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",
                      colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),sep="") -> Title0
              }else
                {
                  paste("Heatmap_",FileNameGeneset,"_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
                  paste("Heatmap_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),sep="") -> Title0
              }
              Title0 %>% paste("\n\n",FileNameGeneset) -> Title1
              file.path(Path0,DirName1,DirNameN0,FileName0) -> FileName1
              # Symbol size
              5000 -> maxheatmap
              c(5,10,20,30,40,50,100,200,1000,maxheatmap) -> NRowheatmap0
              c(1,0.9,0.8,0.7,0.6,0.5,0.35,0.25,0.2,0.1) -> CexRow0
              if(nrow(ttt) < maxheatmap)
              {nrow(ttt) %>% "<="(NRowheatmap0) %>% which %>% min %>% CexRow0[.] -> CexRow1}else{ CexRow0[10] -> CexRow1 }
              pdf(FileName1)
              ttt %>% heatmap(Factor0,labCol=colnames(ttt),labRow=tt$Symbol,
                              dendrogram="row",cexRow=CexRow1,cexCol=0.5)
              title(main=Title1,cex.main=0.50)
              graphics.off()
            }
          }
      }
    if(dopar){parallel::stopCluster(cl); doParallel::stopImplicitCluster()}
    #
  }
  # ----
  # 10 - GMT heatmap all groups
  # ----
  if(!is.null(keywords) & dogmtgenesetheatmap & length(levels(factor)) > 2 & length(t6) > 1)
  {
    length(levels(factor)) -> NbLevels0
    paste("HeatmapGMT_",NbLevels0,sep = "") -> DirName1
    file.path(Path0,DirName1) %>% dir.create
    # moal:::thresholdlist[threshold] -> Threshold0
    # keywords %>% basename %>% sub("(.*).gmt","\\1",.) -> DirNameGMT0
    # file.path(Path0,DirName1,DirNameGMT0) %>% lapply(dir.create)
    # GenesetLists1 <- foreach(i=1:length(keywords), .combine = "c") %do%
    #   {
    #     keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
    #     keywords[i] %>% file("r") -> f
    #     readLines(f,warn=FALSE) -> rl0
    #     close(f)
    #     rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2)) -> GenesetLists0
    #     rl0 %>% strsplit("\t") %>% lapply("[",1) %>% unlist -> Names0
    #     paste(Names0,"|",Keyword0,sep="") -> Names1
    #     names(GenesetLists0) <- Names1
    #   }

    keywords %>% basename %>% sub("(.*).gmt","\\1",.) -> DirNameGMT0
    file.path(Path0,DirName1,DirNameGMT0) %>% lapply(dir.create)
    
    GenesetLists1 <- foreach(i=1:length(keywords), .combine = "c") %do%
      {
        keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
        keywords[i] %>% file("r") -> f
        readLines(f,warn=FALSE) -> rl0
        close(f)
        rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2)) -> GenesetLists0
        rl0 %>% strsplit("\t") %>% lapply("[",1) %>% unlist -> Names0
        paste(Names0,"|",Keyword0,sep="") -> Names1
        names(GenesetLists0) <- Names1
      }
    GenesetLists1 %>% length
    GenesetListsSymbol <- foreach(i=1:length(keywords), .combine = "c") %do%
      {
        keywords[i] %>% basename %>% sub("(.*).gmt","\\1",.) -> Keyword0
        keywords[i] %>% file("r") -> f
        readLines(f,warn=FALSE) -> rl0
        close(f)
        rl0 %>% strsplit("\t") %>% lapply("[",-c(1,2))
      }
    GenesetListsSymbol %>% head
    GenesetListsSymbol %>% length
    GenesetListsSymbol %>% lapply(length)
    
    
    
    
    
    
    
        #
    MatHeatmap0 <- foreach(j=1:length(t6), .combine = "rbind") %do%
      {
        # threshold
        if(THRESHOLD)
        {
          t6[[j]] %>% dplyr::arrange(.data[[colnames(t6[[j]])[2]]]) %>%
            dplyr::filter((.data[[colnames(t6[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                             .data[[colnames(t6[[j]])[3]]] > unlist(Threshold0[[1]][2]) )|(
                               .data[[colnames(t6[[j]])[2]]] < unlist(Threshold0[[1]][1]) & 
                                 .data[[colnames(t5[[j]])[3]]] < -unlist(Threshold0[[1]][2]))) -> p
          
        }else{ t6[[j]] -> p }
        if(species=="hs"){p %>% dplyr::select(1,4) %>% inner_join(dat) -> t }else{p %>% dplyr::select(1,6) %>% inner_join(dat) -> t}
      }
    MatHeatmap0 %>% unique -> MatHeatmap1
    if(dopar){ parallel::detectCores() -> nb ; parallel::makeCluster(nb) -> cl; doParallel::registerDoParallel(cl)}
    foreach(i=1:length(GenesetLists1),.packages=c("magrittr","dplyr","moal","foreach","stringr","igraph","Rgraphviz","moalstringdbhs")) %do%
      {
        # GenesetLists1[[i]] -> symbollist
        GenesetListsSymbol[[i]] -> symbollist
        symbollist %>% annot(.,idtype = "SYMBOL") -> a0
        a0$Symbol -> nodelist
        
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameGeneset
        GenesetLists1[i] %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        
        
        
        # GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> title
        # GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(1) -> FileNameGeneset
        # GenesetLists1[i] %>% names %>% strsplit("\\|") %>% unlist %>% "["(2) -> DirNameN0
        nodelist %>% data.frame(Symbol=.) %>% inner_join(MatHeatmap1) -> tt
        if(nrow(tt) > 2)
        {
          tt %>% dplyr::select(-c(1,2)) -> ttt
          factor -> Factor0
          if(THRESHOLD)
          {
            paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",FileNameGeneset,"_",
                  colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(ttt),".pdf",sep="") -> FileName0
            paste("Heatmap_p",unlist(Threshold0[[1]][3]),"_fc",unlist(Threshold0[[1]][4]),"_",
                  colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(ttt),sep="") -> Title0
          }else{
            paste("Heatmap_",FileNameGeneset,"_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(tt),".pdf",sep="") -> FileName0
            paste("Heatmap_",colnames(t6[[j]])[2] %>% gsub("p_","",.),"_",nrow(ttt),sep="") -> Title0
          }
          Title0 %>% paste("\n\n",FileNameGeneset) -> Title1
          file.path(Path0,DirName1,DirNameN0,FileName0) -> FileName1
          pdf(FileName1)
          ttt %>% heatmap(Factor0,labCol=colnames(ttt),labRow=tt$Symbol,
                          dendrogram="row",cexRow=0.7,cexCol=0.5)
          title(main=Title1,cex.main=0.50)
          graphics.off()
        }    
      }
    if(dopar){parallel::stopCluster(cl); doParallel::stopImplicitCluster()}
    #
  }
  paste("ena done\n",sep="") %>% cat
}
