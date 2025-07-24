MOAL: Multi-Omic Analysis at Lab. A simplified workflow function to make reproducible omic bioanalysis.

MOAL install (v 1.2.1):

install.packages("moal",repos=c("https://fdumbioinfo.r-universe.dev","https://cloud.r-project.org","https://bioc.r-universe.dev"))

annotation depends packages:

install.packages('moalannotgene',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalannotensg',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalannotenst',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalannotensp',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalstringdbhs',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalstringdbmm',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalstringdbrn',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalstringdbdr',repos=c('https://fdumbioinfo.r-universe.dev'))
install.packages('moalstringdbss',repos=c('https://fdumbioinfo.r-universe.dev'))

Workflow summary:
 - Quality controls: histogram, boxplot, PCA and sample hierarchical clustering.
 - Supervised analysis: analysis of variance (ANOVA), Fratio barplot and filtering.
 - Unsupervised analysis for selected features: row hierarchical clustering, PCA and pattern search across factor levels.
 - Graph generation for selected features: volcanoplots, heatmaps, lineplots, boxplots, PCA, Fratio.
 - Functional analysis: GSEA MSigDB enrichment analysis and StringDB interaction network


RopenSci r-universe:
https://fdumbioinfo.r-universe.dev/moal

bioRxiv preprint:
https://doi.org/10.1101/2023.10.17.562686

Zenodo repository:
https://zenodo.org/records/15309968

github:
https://github.com/fdumbioinfo/moal

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

