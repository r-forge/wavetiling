\name{plot}
\alias{plot}

\title{plot model fit and genomic regions}

\description{Plot function to visualize the results of the wavelet-based transcriptome analysis. Both the model fit and the significant genomic regions can be plotted and compared with the annotation.}

\usage{plot(fit, inf, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{\dots}{extra arguments; see \code{\link{plotWfm}}}
}

\details{The plot utilities of the \code{GenomeGraphs}-package constitute the backbone of the \code{plot} function.}

\value{nothing returned}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\seealso{\code{\link{plotWfm}}}

\examples{
\dontrun{
  library(waveTilingData)
  data(tair9gff)
  data(leafdevFit)
  data(leafdevInfCompare)
  gene1 <- tair9gff[tair9gff$ID %in% "AT1G04350",]
  start <- gene1$start-2500
  end <- gene1$end+2500
  plot(fit=leafdevFit,inf=leafdevInfCompare,annoFile=gene1,minPos=start,maxPos=end,two.strand=TRUE,plotData=TRUE,plotMean=FALSE,tracks=c(1,2,6,10,11))
}
}

\keyword{hplot}
