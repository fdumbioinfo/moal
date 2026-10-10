<b><h3>MOAL: Multi Omic Analysis at Lab</h3></b>

Workflow summary:
 - Quality controls: histogram, boxplot, PCA, HC
 - Supervised analysis: analysis of variance
 - volcanoplots, heatmaps, lineplots, boxplots (Kruskal-Wallis), ANOVA Fratio
 - Functional analysis: GSEA MSigDB enrichment analysis, StringDB interaction networks



<b>Install moal and test omic() function workflow:</b>

```r
# -----
# Paris-Saclay - GPL3 - https://github.com/fdumbioinfo/moal
# title: moal omic() demo
# date: 10-10-2026
# -----
#
# Install moal
if(!require("moal",quietly=TRUE)){source("https://raw.githubusercontent.com/fdumbioinfo/rtools/main/moal-demo/0-moal-install-r-universe.r")}
if(!require("data.table",quietly=TRUE)){install.packages("data.table")}
# Test omic() 
library(moal)
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
dat$rowID %>% annotdb(species = "hs", idtype = "GENE", dboutput = "ncbi") -> annot
annot %>% head
annot %>% str
annot %>% dplyr::select(GeneID,Symbol) -> annot
annot %>% head
# omic analysis with MSigDB GSEA functional analysis
dat %>% data.frame -> dat
sif %>% data.frame -> sif
"1-omic-outputdata" -> path
"GSE65055" -> dirname
NULL -> sample
#
omic(
  #
  # input
  dat = dat , sif = sif,
  annot = annot, species = "hs",
  #
  # anova model
  model = "ANEUPLOIDY", padj = "none", anovafc = "anovafc",
  batch = "TISSUE", 
  # paired = NULL,
  #
  # anova all filtering
  # threshold = c(1, 2, 3, 4, 9, 10, 11, 12),
  dovenn = T, docluster = T, dopattern = T,
  #
  # plot
  doheatmap = T, dovolcanoplot = T, dolineplot = T, doboxplotrow = T,
  #
  # functional analysis doena for GSEA and doenaora for ORA significant feature lists
  doena = T, doenaora = F, gmtfiles = NULL,
  dotopheatmap = T, dotopgenesetheatmap = T, dogmtgenesetheatmap = F, 
  dotopnetwork = T, dotopgenesetnetwork = F, dogmtgenesetnetwork = F, 
  # 
  # other par
  path = path, dirname = dirname, sample = sample,
  # parallelisation
  dopar = T, nbcore = NULL, liblocpath = NULL, savehistory = T )
#
```


<p><b>Test <a href="https://github.com/fdumbioinfo/rtools/tree/main/moal-demo">other examples </a> for functional analysis, volcanoplot, heatmap...</b></p>


UMS IPSIT BIOINFO at Paris Saclay University: https://www.ipsit.universite-paris-saclay.fr/?-bioinfo-

r-universe: https://fdumbioinfo.r-universe.dev/moal

bioRxiv preprint: https://doi.org/10.1101/2023.10.17.562686

Zenodo repository: https://zenodo.org/records/15309968

Release updates:

1.2.3

Add replacena function to replace NA or 0 value matrix by background noise distribution

Add normalization method to norm function : TMM, TMMwsp, RLE, upperquartile from EdgeR

Test omic function witout data to simulate random analysis

Change annot function name to annotdb

Automatically annot when dat matrix first column contains Ensembl gene or NCBI gene ID

Automatically join data and annotation tables (so could have different row numbers)

1.2.2

Automatically annotate and make functional analysis when Ensembl or NCBI gene ID provided in first column of dat argument

annotation and functional analysis available for Homo sapiens, Mus musculus, Rattus norvegicus, Sus scrofa, Pan troglodytes, Bos taurus, Drosophila melanogaster, Caenorhabditis elegans

annot argument data.frame can have different number of row compare to dat data.frame

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

