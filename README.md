<b><h3>MOAL: Multi Omic Analysis at Lab</h3></b>

Workflow summary:
 - Quality controls: histogram, boxplot, PCA, HC
 - Supervised analysis: analysis of variance
 - volcanoplots, heatmaps, lineplots, boxplots (Kruskal-Wallis), ANOVA Fratio
 - Functional analysis: GSEA MSigDB enrichment analysis, StringDB interaction networks


<b>Install and test moal::omic() workflow:</b>

```r
# -----
# UMS IPSIT BIOINFO - Licence GPL-3
# https://github.com/fdumbioinfo/moal
# title: moal omic() demo
# date: 11-12-2025
# -----
#
# libraries
if(!require("moal",quietly=TRUE)){source("https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/0-moal-install-r-universe.r")}
if(!require("data.table",quietly=TRUE)){install.packages("data.table")}
library(moal);moal::env()
??moal::omic
# output directory
if(!file.exists("1-omic-outputdata")){"1-omic-outputdata" %>% dir.create}
# loading data
# normalized expression data
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/GSE65055_datanorm_22_23786.tsv" -> url
data.table::fread(url) %>% data.frame -> dat
dat %>% head
# loading sample information file
"https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/inputdata/GSE65055_metadata_4_22.tsv" -> url
data.table::fread(url) %>% data.frame -> sif
sif %>% head
# Ordering factors for pairwise comparison fold-changes
sif$ANEUPLOIDY %>% ordered(c("Control","T13","T18","T21")) -> sif$ANEUPLOIDY
sif$TISSUE %>% as.factor -> sif$TISSUE
# create annotation
dat$rowID %>% moal::annot(species ="hs", idtype = "GENE", dboutput = "ncbi") -> annot
annot %>% head
annot %>% str
annot %>% dplyr::select(GeneID,Symbol) -> annot
annot %>% head
# omic analysis with MSigDB GSEA functional analysis
dat %>% data.frame -> dat
sif %>% data.frame -> sif
moal::omic(dat=dat,sif=sif,annot=annot,species="hs",model="ANEUPLOIDY",batch="TISSUE",dirname="GSE65055",path="1-omic-outputdata")
#
```

<p><b>Test <a href="https://github.com/fdumbioinfo/rtools/tree/main/moal-demo">other examples </a> for functional analysis, volcanoplot, heatmap...</b></p>


UMS IPSIT BIOINFO at Paris Saclay University: https://www.ipsit.universite-paris-saclay.fr/?-bioinfo-

r-universe: https://fdumbioinfo.r-universe.dev/moal

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

