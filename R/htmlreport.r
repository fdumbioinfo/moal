#' @title html report 
#' @description create html report
#' @param omicfolder character path
#' @param name character name report
#' @param dpihtml numerical plot resolution for html
#' @return report folder
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom quarto quarto_render
#' @importFrom htmltools htmlEscape
#' @importFrom pdftools pdf_convert
#' @importFrom foreach foreach
htmlreport <- function( omicfolder = NULL, name = NULL , dpihtml = 150 )
{
  i=1
  omicfolder %>% basename -> dossier_omic
  omicfolder -> Path0
  name -> nomTest
  #
  # create lightbox extension
  # Path0 %>% file.path("_extensions") %>% dir.create
  # Path0 %>% file.path("_extensions","quarto-ext") %>% dir.create
  # Path0 %>% file.path("_extensions","quarto-ext","lightbox") %>% dir.create
  # Path0 %>% file.path("_extensions","quarto-ext","lightbox","resources") %>% dir.create
  # Path0 %>% file.path("_extensions","quarto-ext","lightbox","resources","css") %>% dir.create
  # Path0 %>% file.path("_extensions","quarto-ext","lightbox","resources","js") %>% dir.create
  #
  # moal:::lbext -> lbext
  # i=1
  # foreach(i=1:length(lbext)) %do%
  #   {
  #     Path0 %>% file.path(lbext[[i]][[1]]) %>% file("w") -> fo 
  #     # rl0 %>% writeLines(fo)
  #     cat(paste(lbext[[i]][[2]],collapse="\n"),file=fo)
  #     fo %>% close
  #   }
  # html output
  Path0 %>% file.path("html") %>% dir.create
  #
  # pdf convert to png
  pdfs <- list.files(Path0,pattern="\\.pdf$",full.names=T,ignore.case=T,recursive=T)
  pdfs %>% dirname %>% unique %>% file.path("png_temp") %>% lapply(dir.create, showWarnings = F)
  #
  foreach(i=1:length(pdfs),.packages = c("magrittr","pdftools")) %dopar%
    {
      pdfs[i] %>% dirname %>% file.path("png_temp",pdfs[i] %>% basename %>% gsub("\\.pdf$",".png",.)) -> filenames
      if(file.exists(filenames)){return(NULL)}else{pdftools::pdf_convert(pdf=pdfs[i],format="png",filenames=filenames,dpi=dpihtml)}
    }
  #
  # input omic folder
  dossier_inputdata <- list.dirs(Path0,full.names=T,recursive=F)
  dossier_inputdata <- dossier_inputdata[grepl("inputdata",dossier_inputdata)]
  fichiers_inputdata <- list.files(pattern="\\.tsv$",dossier_inputdata,full.names=T)
  fichier_metadata <- fichiers_inputdata[grepl("^metadata", basename(fichiers_inputdata))]
  dossier_QCsansRBE <- list.dirs(Path0,full.names=T,recursive=F)
  dossier_QCsansRBE <- dossier_QCsansRBE[grepl("^QC_all_",basename(dossier_QCsansRBE)) & !grepl("RBE",basename(dossier_QCsansRBE))][1]
  dossier_QCRBE <- list.dirs(Path0,full.names=T,recursive=F)
  dossier_QCRBE <- dossier_QCRBE[grepl("^QC_all",basename(dossier_QCRBE)) & grepl("RBE",basename(dossier_QCRBE))][1]
  fichier_HistosansRBE <- list.files(file.path((dossier_QCsansRBE),"png_temp"))
  fichier_HistoRBE <- list.files(file.path((dossier_QCRBE),"png_temp"))
  dossier_anovastats <- Filter(function(x) grepl("anovastats",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  fichier_anovastats_all <- list.files(pattern="\\.pdf$",dossier_anovastats)
  fichier_anovastats_all <- fichier_anovastats_all[grepl("anovastats_all",basename(fichier_anovastats_all))]
  dossier_venn <- Filter(function(x) grepl("venn",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_DiffLists <- Filter(function(x) grepl("DiffLists_", basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_patternpval <- Filter(function(x) grepl("patternpval_",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_patternfc <- Filter(function(x) grepl("patternfc_",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_cluster <- Filter(function(x) grepl("cluster",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_volcanoplot <- Filter(function(x) grepl("volcanoplot",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_pca <- Filter(function(x) grepl("pca",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_lineplot <- Filter(function(x) grepl("lineplot_",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_boxplot <- Filter(function(x) grepl("boxplot_",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  dossier_ena <- Filter(function(x) grepl("ena_",basename(x)),list.dirs(Path0,full.names=T,recursive=F))
  #
  # qmd script to render
  contenu <- r"---(---
title: "__NOMTEST__"

params:
  famille: "tout"
  dossier: ""

execute:
  echo: false
  warning: false
  message: false
  error: true

toc: __TOC__
toc-depth: 4

format:
  html:
  #embarquer les images dans le fichier HTML. Pour pouvoir effacer les dossiers avec les pdf tansformes en png.
    lightbox: auto
    embed-resources: true
    number-sections: true
    number-depth: 4
---
  
  
```{=html}
<style>
  
.methode-3col {
  display: flex;
  flex-wrap: wrap;  
  justify-content: flex-start;
  align-items: flex-start;
  gap: 1%;
}
 
.methode-3col > div {
  width: 30%;   
  box-sizing: border-box;
}
.methode-3col > div img {
  width: 90%; /*image occupe 90% de la colonne*/
  height: auto;  
  display: block;
  margin: 0 auto;   /* centre l'image dans sa colonne */
}
.methode-3col figure {
  margin: 0;
  width: 90%;          /* meme largeur que l image */
  margin-left: auto;
  margin-right: auto;
  display: flex;
  flex-direction: column;
  align-items: center;
}
/*enleve les marges superflues ici entre mes images*/
.methode-3col figcaption {
  font-size: clamp(8px, 0.9vw, 11px);   /* adapte a une colonne etroite */
    line-height: 1.25;                  /* evite le chevauchement entre lignes de texte */
    text-align: center;
  word-wrap: break-word;
  overflow-wrap: break-word;
  white-space: normal;
  max-width: 100%;
  margin: 4px auto 0;
  padding: 0 4px;
}
 
figcaption, .figure-caption {
  text-align: center;
  font-size: clamp(10px, 1vw, 13px);
  word-wrap: break-word;
  overflow-wrap: break-word;
  white-space: normal;
  max-width: 100%;
  margin: 0 auto;
  padding: 0 4px;
}
 
h1, h2, h3, h4, h5, h6, .title {
  overflow-wrap: anywhere;
  word-break: break-word;
  max-width: 100%;
  box-sizing: border-box;
}

h2 { padding-left: 1rem; }
h3 { padding-left: 2rem; }
h4 { padding-left: 3rem; }
h5 { padding-left: 4rem; }
h6 { padding-left: 5rem; }
/* Le contenu (element de grille) peut retrecir au lieu d etre pousse par un titre long */
#quarto-content, main.content, .page-columns .content {
  min-width: 0;
}

#TOC .header-section-number,
#TOC .toc-section-number,
nav[role="doc-toc"] .header-section-number,
nav[role="doc-toc"] .toc-section-number,
#quarto-margin-sidebar .header-section-number,
#quarto-margin-sidebar .toc-section-number {
  display: none;
}

#TOC, #quarto-margin-sidebar, nav[role="doc-toc"] {
  overflow-x: hidden;
}
#TOC a, nav[role="doc-toc"] a, #quarto-margin-sidebar .nav-link {
  white-space: normal;
  overflow-wrap: anywhere;
  word-break: break-word;
}
  
pre {
  white-space: pre-wrap;       /* retour a la ligne automatique */
  word-wrap: break-word;       /* coupe les mots trop longs si besoin */
  overflow-x: hidden;          /* supprime la barre de defilement horizontale */
}
pre, pre code {
  white-space: pre-wrap !important;
  word-break: break-all !important;
  overflow-wrap: break-word !important;
  overflow-x: hidden !important;
}
</style>
```

```{r setup}
#| include: false

library(magrittr)
library(pdftools)
library(magick)
library(gtools)
library(foreach)
# library(parallel)
# library(doParallel)
# library(quarto)
```

```{r}
#| label: chemins
#| include: false

nomTest <- "__NOMTEST__"
dossier_omic <- __DOSSIER_OMIC__
# nbcore <- __NBCORE__

dossier_inputdata <- list.dirs(dossier_omic,full.names=T,recursive=F)
dossier_inputdata<- dossier_inputdata[grepl("inputdata",dossier_inputdata)]
fichiers_inputdata <- list.files (pattern = "\\.tsv$", dossier_inputdata, full.names = T)
fichier_metadata <- fichiers_inputdata[grepl("^metadata", basename(fichiers_inputdata))]
dossier_QCsansRBE <- list.dirs(dossier_omic,full.names=T,recursive=F)
dossier_QCsansRBE <- dossier_QCsansRBE[grepl("^QC_all_",basename(dossier_QCsansRBE)) & !grepl("RBE",basename(dossier_QCsansRBE))][1]
dossier_QCRBE <- list.dirs(dossier_omic,full.names=T,recursive=F)
dossier_QCRBE <-dossier_QCRBE[grepl("^QC_all",basename(dossier_QCRBE)) & grepl("RBE",basename(dossier_QCRBE))][1]
fichier_HistosansRBE <- list.files(file.path((dossier_QCsansRBE),"png_temp"))
fichier_HistoRBE <- list.files(file.path((dossier_QCRBE),"png_temp"))
dossier_anovastats <- Filter(function(x) grepl("anovastats",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
fichier_anovastats_all <- list.files(pattern="\\.pdf$",dossier_anovastats)
fichier_anovastats_all <- fichier_anovastats_all[grepl("anovastats_all",basename(fichier_anovastats_all))]
dossier_venn <- Filter(function(x) grepl("venn",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_DiffLists <- Filter(function(x) grepl("DiffLists_", basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_patternpval <- Filter(function(x) grepl("patternpval_",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_patternfc <- Filter(function(x) grepl("patternfc_",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_cluster <- Filter(function(x) grepl("cluster",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_volcanoplot <- Filter(function(x) grepl("volcanoplot",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_pca <- Filter(function(x) grepl("pca",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_lineplot <- Filter(function(x) grepl("lineplot_",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_boxplot <- Filter(function(x) grepl("boxplot_",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
dossier_ena <- Filter(function(x) grepl("ena_",basename(x)),list.dirs(dossier_omic,full.names=T,recursive=F))
#
```

```{r}
#| label: affichage-lien
#| results: 'asis'
#| echo: false

afficher_lien <- function(dossier){
  if (is.null(dossier) || length(dossier) == 0 || is.na(dossier)) return(invisible(NULL))
  nom_lien <- basename(dossier)
  chemin <- gsub("\\\\", "/", normalizePath(dossier))
  chemin <- utils::URLencode(sub("^/", "", chemin))
  cat(sprintf("[%s](file:///%s)", nom_lien, chemin))
}
```

```{r}
#| echo: false
#| label: affichage-index
#| results: 'asis'

familles_dossiers <- list(
  inputdata   = dossier_inputdata,
  QC          = Filter(Negate(is.na), c(dossier_QCsansRBE, dossier_QCRBE)),
  anovastats  = dossier_anovastats,
  venn        = dossier_venn,
  DiffLists   = dossier_DiffLists,
  patternpval = dossier_patternpval,
  patternfc   = dossier_patternfc,
  cluster     = dossier_cluster,
  volcanoplot = dossier_volcanoplot,
  pca         = dossier_pca,
  lineplot    = dossier_lineplot,
  boxplot     = dossier_boxplot,
  ena         = dossier_ena
)

familles_dossiers <- Filter(function(x) length(x) > 0, familles_dossiers)

ordre_dossiers <- if (nzchar(params$dossier)) {
  params$dossier
} else if (params$famille == "tout") {
  unlist(familles_dossiers, use.names = FALSE)
} else {
  familles_dossiers[[params$famille]]
}

ordre_dossiers <- ordre_dossiers[dir.exists(ordre_dossiers)]

dossiers_ordonnes_afficher <- function(dossiers){
  for (d in dossiers) {
    dossier_afficher(d, 1)
  }
}
```

```{r}
#| label: lien-retour
#| results: 'asis'

cat(sprintf('<p><a href="%s"> return to index</a></p>\n',paste0("../",nomTest,"_index.html")))
cat("<br>\n\n")
```

```{r}
#| label: affichage-images-ss-dossiers
#| results: 'asis'

# strrep("#", niveau) construit le bon nombre de # (niveau 1 #, niveau 2 ##, etc.).

dossier_afficher<- function(dossier, niveau){
  
  cat("\n\n", strrep("#", min(niveau, 6)), " ", basename(dossier), "\n\n", sep = "")
  
  afficher_lien(dossier)
  cat("<br>\n\n")
  
  tous_fichiers <- list.files(file.path(dossier, "png_temp"), pattern = "\\.png$", full.names = T,  recursive = F)

  if (length(tous_fichiers) == 1 && grepl("histo", basename(tous_fichiers), ignore.case = T)) {
    cat(sprintf('\n![%s](%s){width=70%% .lightbox}\n\n',
              tools::file_path_sans_ext(basename(tous_fichiers)),
              tous_fichiers))
    } else if (length(tous_fichiers) > 0) {
    afficher3col(dossier)
  }
  
  sous_dossiers <- list.dirs(dossier, full.names = T, recursive = F)
  
  sous_dossiers <- sous_dossiers[basename(sous_dossiers) != "png_temp"]
  
  if (length(sous_dossiers) != 0) {
    ordre <- gtools::mixedorder(basename(sous_dossiers))
    dossiers_ordonnes <- sous_dossiers[ordre]
    for (d in dossiers_ordonnes) {
      
      dossier_afficher(d, niveau = niveau + 1)
      
    }
  }
}
```

```{r}
#| echo: false
#| results: 'asis'
#| label: afficher3col

afficher3col <- function(dossier){
  
  tous_fichiers <- list.files(file.path(dossier, "png_temp"), pattern = "\\.png$", full.names = T)  
  ordre <- gtools::mixedorder(basename(tous_fichiers))
  fichiers_ordonnes <- tous_fichiers[ordre]
  ordre2 <- order(!grepl("Venn", fichiers_ordonnes))
  fichiers_ordonnes2 <- fichiers_ordonnes[ordre2]
  nomf <- tools::file_path_sans_ext(basename(fichiers_ordonnes2))
  
  cat("\n")
  
  classe <-".methode-3col"
  cat(sprintf("\n::: {%s}\n\n", classe))
  cat(sprintf("![%s](%s){.lightbox}\n\n",
              nomf,
              fichiers_ordonnes2), sep = "")
  cat(":::\n")
}

```

```{r}
#| label: Afficher-log

fichier_log <- file.path(dossier_omic, "log.txt")
if (params$famille %in% c("tout", "inputdata") && file.exists(fichier_log)) {
  cat(readLines(fichier_log, warn = FALSE), sep = "\n")
}
```


```{r}
#| label: afficher-tout
#| results: 'asis'

if (params$famille %in% c("tout", "inputdata")) {
  afficher_lien(dossier_omic)
}
dossiers_ordonnes_afficher(ordre_dossiers)
```

)---"
 #
 # -----
 # 1 - omic folder path for qmd and report name
 # -----
 #
 contenu <- gsub("__NOMTEST__",nomTest,contenu,fixed=T)
 contenu <- gsub("__DOSSIER_OMIC__",deparse("./"),contenu,fixed=T)
 # if(!is.null(nbcore)){ contenu <- gsub("__NBCORE__",nbcore,contenu,fixed=T) }else{ contenu <- gsub("__NBCORE__","NULL",contenu,fixed=T) }
 #
 # -----
 # 2 - create qmd
 # -----
 #
 for(avec_toc in c(TRUE, FALSE))
 {
   writeLines(
     enc2utf8(gsub("__TOC__",tolower(avec_toc),contenu,fixed=T)),
     paste0(nomTest, if (avec_toc) "_toc" else "_notoc",".qmd") %>% file.path(Path0,.),useBytes=TRUE)
 }
 # -----
 # 3 - render
 # -----
 #
 # one html by folder
 #
 # sous_dossiers_de <- function(dossiers)
 # {
 #   sd <- list.dirs(dossiers, full.names = T, recursive = F)
 #   sd <- sd[basename(sd) != "png_temp"]
 #   if (length(sd) > 1) sd else dossiers
 # }
 #
 libelles <- c(
   inputdata   = "inputdata",
   QC          = "Quality Controls (QC)",
   anovastats  = "ANOVA",
   venn        = "Diagrammes de Venn",
   DiffLists   = "DiffLists",
   patternpval = "Patterns p-value",
   patternfc   = "Patterns fc",
   cluster     = "Cluster",
   volcanoplot = "Volcanoplot",
   pca         = "PCA",
   lineplot    = "Lineplot",
   boxplot     = "Boxplot",
   ena         = "ENA"
 )
 #
 familles_dossiers <- list(
   inputdata   = dossier_inputdata,
   QC          = Filter(Negate(is.na), c(dossier_QCsansRBE, dossier_QCRBE)),
   anovastats  = dossier_anovastats,
   venn        = dossier_venn,
   DiffLists   = dossier_DiffLists,
   patternpval = dossier_patternpval,
   patternfc   = dossier_patternfc,
   cluster     = dossier_cluster,
   volcanoplot = dossier_volcanoplot,
   pca         = dossier_pca,
   # lineplot    = sous_dossiers_de(dossier_lineplot),
   # boxplot     = sous_dossiers_de(dossier_boxplot),
   lineplot    = dossier_lineplot,
   boxplot     = dossier_boxplot,
   ena         = dossier_ena )
 #
 familles_dossiers <- Filter(function(x) length(x) > 0, familles_dossiers)
 rapports <- do.call(rbind, lapply(seq_along(familles_dossiers), function(i) {
  chemins <- unique(as.character(unlist(familles_dossiers[[i]])))
  chemins <- chemins[!is.na(chemins) & dir.exists(chemins)]
  if(length(chemins) == 0) return(NULL)
  data.frame(famille = names(familles_dossiers)[i],
             chemin  = chemins,
             stringsAsFactors = FALSE)
 }))
 #
 if(is.null(rapports) || nrow(rapports) == 0){stop("no folder found : ",Path0)}
 rapports$fichier <- paste0(nomTest,"_",rapports$famille,"_",basename(rapports$chemin),".html")
 nb_titres <- function(dossier)
 {
   d <- list.dirs(dossier, full.names = TRUE, recursive = TRUE)
   sum(!basename(d) %in% c("png_temp"))
 }
 #
 # for(k in seq_len(nrow(rapports))) 
 k=1
 foreach(k=1:nrow(rapports),.packages=(c("quarto"))) %do% 
   {
     paste0(nomTest, if (nb_titres(rapports$chemin[k]) > 1) "_toc" else "_notoc", ".qmd") -> entree
     entree %>% file.path(Path0,.) -> entree
     quarto::quarto_render(
       input = entree,
       execute_params = list(famille = rapports$famille[k],dossier = rapports$chemin[k] %>% basename),
       output_file = rapports$fichier[k] )
    }
 # -----
 # 4 - create index.html
 # -----
 nom_index <- file.path(".",paste0(nomTest,"_index.html"))
 blocshtmlindex <- function(f)
 {
   rapports[rapports$famille == f, ] -> sous
   unname(libelles[f]) -> titre
   if(is.na(titre)) titre <- f
   c(sprintf('<h2>%s</h2>', htmltools::htmlEscape(titre)),'<ul>',
     sprintf('  <li><a href="%s">%s</a></li>',file.path("html",sous$fichier), htmltools::htmlEscape(basename(sous$chemin))),'</ul>')
  }
 rapports$famille %>% unique %>% lapply(blocshtmlindex) %>% unlist -> blocs
 # blocs <- unlist(lapply(unique(rapports$famille), function(f) {
 #   sous  <- rapports[rapports$famille == f, ]
 #   titre <- unname(libelles[f])
 #   if(is.na(titre)) titre <- f
 #   c(sprintf('<h2>%s</h2>', htmltools::htmlEscape(titre)),'<ul>',
 #   sprintf('  <li><a href="%s">%s</a></li>',sous$fichier, htmltools::htmlEscape(basename(sous$chemin))),'</ul>')
 # }))
 #
 html <- c(
   '<!DOCTYPE html>',
   '<html lang="fr"><head><meta charset="utf-8">',
   sprintf('<title>Index %s</title>', htmltools::htmlEscape(basename(dossier_omic))),
   '<style>',
   '  body { font-family: sans-serif; max-width: 70%; margin: 40px auto; padding: 0 16px; }',
   '  h2 { margin: 28px 0 4px; font-size: 1.3em; }',
   '  li { margin: 6px 0; font-size: 1.1em; }',
   '</style></head><body>',
   sprintf('<h1>%s</h1>', htmltools::htmlEscape(basename(dossier_omic))),
   sprintf('<p>create at %s</p>', format(Sys.time(), "%d/%m/%Y %H:%M")),
   blocs,
   '</body></html>' )
 #
 writeLines(html,nom_index %>% file.path(Path0,.) ,useBytes=T)
 # -----
 # 5 - remove temp files
 # -----
 #
 # remove png
 dossiers <- list.dirs(file.path(Path0),recursive=T,full.names=T)
 dossiers_png_temp <- dossiers[basename(dossiers) == "png_temp"]
 unlink(dossiers_png_temp,recursive=T)
 #
 # remove move to html folder
 Path0 %>% list.files -> lhtml0
 lhtml0 %>% grep(".html$",.) %>% lhtml0[.] -> lhtml1
 lhtml1 %>% grep("index.html$",.) -> selindexhtml
 lhtml1[-selindexhtml] %>% file.path(Path0,.) %>% file.copy(Path0 %>% file.path("html"))
 lhtml1[-selindexhtml] %>% file.path(Path0,.) %>% file.remove
 #
 # remove lightbox
 # Path0 %>% file.path("_extensions") %>% list.files(full.names=T,recursive=T) -> llb
 # llb %>% file.remove
 # Path0 %>% file.path("_extensions") %>% unlink(recursive=T)
 #
 # remove qmd
 Path0 %>% list.files -> lqmd0
 lqmd0 %>% grep(".qmd$",.) %>% lqmd0[.] %>% file.path(Path0,.) %>% file.remove
}