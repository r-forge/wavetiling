\name{plotWfm}
\alias{plotWfm}

\title{plot model fit and genomic regions}

\description{Plot function to visualize the results of the wavelet-based transcriptome analysis. Both the model fit and the significant genomic regions can be plotted and compared with the annotation.}

\usage{plotWfm(fit, inf, biomartObj, minPos, maxPos, trackFeature="exon", two.strand=TRUE, plotData=TRUE, plotMean=TRUE, tracks=0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{biomartObj}{object of class TranscriptDb representing an annotation database generated from BioMart.}
  \item{minPos}{minimum genomic position to plot}
  \item{maxPos}{maximum genomic position to plot}
  \item{trackFeature}{track feature. See \code{GenomeGraphs}-package. Default is \code{exon}.}
  \item{two.strand}{logical indicating whether to plot two strands or not}
  \item{plotData}{logical indicating whether to plot the raw data or not}
  \item{plotMean}{logical indicating whether to plot the fitted overall mean function or not}
  \item{tracks}{vector of integers containing the track numbers to plot. Track numbers correspond with the order of the elements in the list output from the \code{getGenmicRegions}-function.}
}

\details{The plot utilities of the \code{GenomeGraphs}-package constitute the backbone of the \code{plotWfm} function.}

\value{nothing returned}

\references{
[1] Durinck S, Bullard J, Spellman PT, Dudoit S: GenomeGraphs: integrated genomic data visualization with R.
BMC Bioinformatics 2009, 10:Article 2.
}


\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  library(waveTilingData)
  library(TxDb.Athaliana.BioMart.plantsmart12)
  data(leafdevFit)
  data(leafdevInfCompare)
  trs <- transcripts(TxDb.Athaliana.BioMart.plantsmart12)
  sel <- trs[elementMetadata(trs)$tx_name \%in\% "AT1G62500.1",]
  start <- start(ranges(sel))-2000
  end <- end(ranges(sel))+2000
  plotWfm(fit=leafdevFit,inf=leafdevInfCompare,biomartObj=TxDb.Athaliana.BioMart.plantsmart12,minPos=start,maxPos=end,two.strand=TRUE,plotData=TRUE,plotMean=FALSE,tracks=c(1,2,6,10,11))
}

\keyword{hplot}
