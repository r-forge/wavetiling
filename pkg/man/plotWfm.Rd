\name{plotWfm}
\alias{plotWfm}

\title{plot model fit and genomic regions}

\description{Plot function to visualize the results of the wavelet-based transcriptome analysis. Both the model fit and the significant genomic regions can be plotted and compared with the annotation.}

\usage{plotWfm(fit, inf, annoFile, minPos, maxPos, trackFeature="exon", overlayFeature=c("gene","transposable_element_gene"), two.strand=TRUE, plotData=TRUE, plotMean=TRUE, tracks=0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{annoFile}{annotation file. This annotation file should be a data.frame with (at least) the following information: \code{chromosome} giving the chromosme of the annotated region; \code{strand} providing the strand orientation of the annotated region; \code{start} and \code{end} containing the start and end position of the annotated region; \code{feature} giving the type of the annotated region (gene, pseudogene, exon,...); \code{ID} containing the annotation ID of the annotated region.}
  \item{minPos}{minimum genomic position to plot}
  \item{maxPos}{maximum genomic position to plot}
  \item{trackFeature}{track feature. See \code{GenomeGraphs}-package. Default is \code{exon}.}
  \item{overlayFeature}{overlay feature. See \code{GenomeGraphs}-package. Default is \code{gene}.}
  \item{two.strand}{logical indicating whether to plot two strands or not}
  \item{plotData}{logical indicating whether to plot the raw data or not}
  \item{plotMean}{logical indicating whether to plot the fitted overall mean function or not}
  \item{tracks}{vector of integers containing the track numbers to plot. Track numbers correspond with the order of the elements in the list output from the \code{getGenmicRegions}-function.}
}

\details{The plot utilities of the \code{GenomeGraphs}-package constitute the backbone of the \code{plotWfm} function.}

\value{nothing returned}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(tair9gff)
  data(leafdevFit)
  data(leafdevInfCompare)
  gene1 <- tair9gff[tair9gff$ID %in% "AT1G04350",]
  start <- gene1$start-2500
  end <- gene1$end+2500
  plotWfm(fit=leafdevFit,inf=leafdevInfCompare,annoFile=gene1,minPos=start,maxPos=end,two.strand=TRUE,plotData=TRUE,plotMean=FALSE,tracks=c(1,2,6,10,11))
}
}

\keyword{hplot}
