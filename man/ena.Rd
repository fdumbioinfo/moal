% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ena.r
\name{ena}
\alias{ena}
\title{Gene set enrichment analysis and interaction network}
\usage{
ena(
  omicdata = NULL,
  gmtfiles = NULL,
  species = "hs",
  dat = NULL,
  factor = NULL,
  filtergeneset = NULL,
  threshold = 1,
  topdeg = 50,
  rangedeg = NULL,
  topena = 50,
  twotailena = TRUE,
  topgeneset = 50,
  intmaxdh = 5000,
  nodesize = 0.6,
  bg = 25000,
  doena = TRUE,
  gsearank = "logfc",
  layout = 1,
  dotopnetwork = TRUE,
  dotopgenesetnetwork = FALSE,
  dogmtgenesetnetwork = FALSE,
  dotopheatmap = TRUE,
  dotopgenesetheatmap = FALSE,
  dogmtgenesetheatmap = TRUE,
  path = NULL,
  dirname = NULL,
  dopar = TRUE
)
}
\arguments{
\item{omicdata}{character data.frame see details}

\item{gmtfiles}{character gmt files list path}

\item{species}{character hs mm rn dr ss}

\item{dat}{data.frame file paths}

\item{factor}{factor factor for heatmap color}

\item{filtergeneset}{character list to filter MSigDB geneset collection}

\item{threshold}{numeric pval 0.05 fc 1.5 by default see details}

\item{topdeg}{numeric top feature to plot on network}

\item{rangedeg}{numeric top DEGs from 1 to topdeg by rangedeg to plot on network}

\item{topena}{numeric top geneset for ena plot}

\item{twotailena}{logical do}

\item{topgeneset}{numeric top geneset number to plot on network}

\item{intmaxdh}{numeric maximum number of interaction to use for Davidson and Harel algorithm layout}

\item{nodesize}{numeric change Symbol size}

\item{bg}{numeric background used for functional analysis over-representation test}

\item{doena}{logical do MSigDB enrichment analysis}

\item{gsearank}{character to choose gsea rank type among fc (by default) logration logfc sqrt}

\item{layout}{numeric for layout neetwork 1 fr by default 2 dh 3 tree 4 circle 5 grid 6 sphere}

\item{dotopnetwork}{logical do top networks}

\item{dotopgenesetnetwork}{logical do geneset networks}

\item{dogmtgenesetnetwork}{logical do keyword networks}

\item{dotopheatmap}{logical do top heatmap}

\item{dotopgenesetheatmap}{logical do geneset heatmap}

\item{dogmtgenesetheatmap}{logical do keyword heatmap}

\item{path}{character for relative path of output directory}

\item{dirname}{character name for output}

\item{dopar}{logical TRUE for parallelization}
}
\value{
file with enrichment analysis results
}
\description{
Gene set enrichment analysis and interaction network
}
\details{
omicdata needs a data.frame with at list 4 column: rowID, (p-values,fold-change) x N and Symbol annotation.

Symbol list are accepted to make ORA enrichment analysis.

To generate heatmap dat and factor parameter are needed. dat accepted complete matrix with rowID for first column.

dat row IDs must match with omicdata row IDs.

Make MSigDB enrichment analysis using GSEA method for non filtering list as input (> 2000)

Make MSigDB Over-Representation enrichment analysis (ORA) using Fisher exact test for list < 2000

Generate STRINGDB interaction network and heatmap for top geneset according to topena par (80 by default)

Only features with p-values < 0.05 et fold-change > 1.1 are displayed on geneset heatmaps (threshold = 1 by default).

See omic function details to display all threshold
}
\examples{
# not run
# gsenain( omicdata , species = "mm")
}
\author{
Florent Dumont \href{mailto:florent.dumont@universite-paris-saclay.fr}{florent.dumont@universite-paris-saclay.fr}
}
