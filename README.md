<b>MOAL: Multi Omic Analysis at Lab</b>

Workflow summary:
 - Quality controls: histogram, boxplot, PCA, HC
 - Supervised analysis: analysis of variance
 - volcanoplots, heatmaps, lineplots, boxplots (Kruskal-Wallis), ANOVA Fratio
 - Functional analysis: GSEA MSigDB enrichment analysis, StringDB interaction networks

Install from <a href="https://fdumbioinfo.r-universe.dev/moal">r-universe</a> (v 1.2.1):

```r
# -----
# UMS IPSIT BIOINFO - Licence GPL-3
# https://github.com/fdumbioinfo/moal
# title: moal install from r-universe
# date: 11122025
# -----
# 
# moal install
options(pkgType = "binary")
# annotation packages
if(!require("moalannotgene",quietly=TRUE)){install.packages("moalannotgene",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensg",quietly=TRUE)){install.packages("moalannotensg",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotenst",quietly=TRUE)){install.packages("moalannotenst",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalannotensp",quietly=TRUE)){install.packages("moalannotensp",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbhs",quietly=TRUE)){install.packages("moalstringdbhs",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbmm",quietly=TRUE)){install.packages("moalstringdbmm",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbrn",quietly=TRUE)){install.packages("moalstringdbrn",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbdr",quietly=TRUE)){install.packages("moalstringdbdr",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
if(!require("moalstringdbss",quietly=TRUE)){install.packages("moalstringdbss",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
# depend packages
if(!require("BiocManager",quietly=TRUE)){install.packages("BiocManager")}
if(!require("Rgraphviz",quietly=TRUE)){BiocManager::install("Rgraphviz",update=F)}
if(!require("limma",quietly=TRUE)){BiocManager::install("limma",update=F)}
if(!require("fgsea",quietly=TRUE)){BiocManager::install("fgsea",update=F)}
# moal package
if(!require("moal",quietly=TRUE)){install.packages("moal",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}
#
```


```r
# -----
# UMS IPSIT BIOINFO - Licence GPL-3
# https://github.com/fdumbioinfo/moal
# title: moal ena() demo
# date: 08-12-2025
# -----
#
# -----
# 1 - functional analysis using Over-Representation analysis (ORA) for Symbol list 
# -----
#
# loading libraries
if(!require("moal",quietly=TRUE)){source("https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/0-moal-install-r-universe.r")}
if(!require("data.table",quietly=TRUE)){install.packages("data.table")}
library(moal);moal::env()
??moal::ena
# output directory
if(!file.exists("2-ena-outputdata")){"2-ena-outputdata" %>% dir.create}
# loading data
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/List_p05_fc15_ANEUPLOIDY_ud_T21vsControl_243.tsv" -> url
data.table::fread(url) %>% data.frame -> omicdata
omicdata %>% str
omicdata %>% dplyr::select(rowID,p_T21vsControl,fc_T21vsControl,Symbol) -> omicdata
omicdata %>% str
#
# example with Symbol list without fold-change and p-value Symbol
#
moal::ena(omicdata=omicdata$Symbol,species="hs",dirname="T21vsControl_243",path="2-ena-outputdata")
#
# example with Symbol list with fold-change and p-value
#
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/GSE65055_metadata_4_22.tsv" -> url
data.table::fread(url) %>% data.frame -> sif
sif$ANEUPLOIDY %>% ordered(c("Control","T13","T18","T21")) -> factor
factor
moal::ena(omicdata=omicdata,factor=factor,species="hs",dirname="fc_T21vsControl_243",path="2-ena-outputdata")
#
# example with Symbol list with fold-change and p-value and expression data for heatmaps
#
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/GSE65055_datanorm_RBE_TISSUE_22_23786.tsv" -> url
data.table::fread(url) %>% data.frame -> dat
dat %>% head
moal::ena(omicdata=omicdata,dat=dat,factor=factor,species="hs",dirname="heatmaps_T21vsControl_243",path="2-ena-outputdata")
#
# example with Symbol list with fold-change and p-value and expression data for heatmaps and
# adding gmtfiles for heatmap on specific genesets
#
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/c2_apoptosis_lung.gmt" %>% download.file(destfile = "c2_apoptosis_lung.gmt")
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/c2_apoptosis_skin.gmt" %>% download.file(destfile = "c2_apoptosis_skin.gmt")
list.files(full.names=T) %>% grep(".gmt$",.,value=T) -> gmtfiles
moal::ena(omicdata=omicdata,dat=dat,factor=factor,gmtfiles=gmtfiles,species="hs",dirname="gmtfiles_T21vsControl_243",path="2-ena-outputdata")
#
# -----
# 2 - functional analysis using GeneSet Enrichment analysis (GSEA)
# -----
#
# loading data
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/GSEA_T21vsCTL_23786.tsv" -> url
data.table::fread(url) %>% data.frame -> omicdata
omicdata %>% str
omicdata %>% dplyr::select(rowID,p_T21vsControl,fc_T21vsControl,Symbol) -> omicdata
omicdata %>% str
#
# example with not filtered differential analysis results with fold-change and p-value and expression data for heatmaps
#
moal::ena(omicdata=omicdata,dat=dat,factor=factor,species="hs",dirname="GSEA_T21vsControl_23786",path="2-ena-outputdata")
#
```

<p>
<b>See more examples at <a href="https://github.com/fdumbioinfo/rtools/tree/main/moal-demo">fdumbioin/rtools/moal-demo/</a> repository.</b>
</p>


UMS IPSIT BIOINFO at Paris Saclay University: https://www.ipsit.universite-paris-saclay.fr/?-bioinfo-

bioRxiv preprint: https://doi.org/10.1101/2023.10.17.562686

Zenodo repository: https://zenodo.org/records/15309968

Release updates:

1.2.1

Fix bug for unbalanced design case using anova model with interaction.
choose gsea rank
Add three factors interaction design


1.2.0

Display improvements on graphics.

Improve omic() function workflow time computing.

Improve omic() output directory browsing and volume.

Enrichment function gsneain():
- GSEA enrichment analysis applied to ANOVA results
- Create heatmaps and interaction networks for significative gene sets
- Heatmap created with gene set GMT files provided
- Biological databases updates: NCBI gene, NCBI orthologs, Ensemble gene, Ensemble transcript, Ensemble proteins, MSigDB gene set collections and StringDB interactions.

