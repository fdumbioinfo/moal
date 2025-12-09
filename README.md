MOAL: Multi-Omic Analysis at Lab. 
A simplified workflow function to make reproducible omic bioanalysis.

Workflow summary:
 - Quality controls: histogram, boxplot, PCA and sample hierarchical clustering.
 - Supervised analysis: analysis of variance (ANOVA), Fratio barplot and filtering.
 - Unsupervised analysis for selected features: row hierarchical clustering, PCA and pattern search across factor levels.
 - Graph generation for selected features: volcanoplots, heatmaps, lineplots, boxplots, PCA, Fratio.
 - Functional analysis: GSEA MSigDB enrichment analysis and StringDB interaction network


Install from r-universe (v 1.2.1):
```r
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
# depends packages
if(!require("BiocManager",quietly=TRUE)){install.packages("BiocManager")}
if(!require("Rgraphviz",quietly=TRUE)){BiocManager::install("Rgraphviz",update=F)}
if(!require("limma",quietly=TRUE)){BiocManager::install("limma",update=F)}
if(!require("fgsea",quietly=TRUE)){BiocManager::install("fgsea",update=F)}
# moal
if(!require("moal",quietly=TRUE)){install.packages("moal",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org"))}```
```

```r
# omic() workflow example using internal GEO dataset GSE65055 (doi: 10.1111/cge.12731):
# See help("omic") to copy and run code.
# loading libraries
library(moal);moal::env()
# loading norm data
moal:::GSE65055normdata -> dat
# loading sample information file
moal:::GSE65055sampledata -> sif
# Ordering factors for pairwise comparison fold-changes
sif$ANEUPLOIDY %>% ordered(c("Control","T13","T18","T21")) -> sif$ANEUPLOIDY
sif$TISSUE %>% as.factor -> sif$TISSUE
# create annotation
dat$rowID %>% moal::annot(species= "hs",idtype="GENE",dboutput = "ncbi") -> annot
annot
# omic analysis
moal::omic(dat,sif,annot,species="hs",model="ANEUPLOIDY",paired="TISSUE",dirname="GSE65055")
#
```





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

