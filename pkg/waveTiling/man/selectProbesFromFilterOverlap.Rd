\name{selectProbesFromFilterOverlap}
\alias{selectProbesFromFilterOverlap}

\title{select probes from \code{MapFilterProbe} object}

\description{Extract the probe indices from a \code{MapFilterProbe} object that map to a region between two specified genomic positions}

\usage{selectProbesFromFilterOverlap(object, chromosome, strand=c("forward","reverse"), minPos=min(getPosition(object)), maxPos=max(getPosition(object)))}

\arguments{
  \item{object}{object of class \code{MapFilterProbe}}
  \item{chromosome}{chromosome}
  \item{strand}{strand}
  \item{minPos}{minimum genomic position}
  \item{maxPos}{maximum genomic position}
}

\value{A list of 2 elements is returned. The first element "selection" gives the probe indices in the filtered \code{MapFilterProbe} object; the second element "selectionInit" gives the probe indices in the original \code{WaveTilingFeatureSet} object.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  library(waveTilingData)
  data(leafdevMapAndFilterTAIR9)
  tt <- selectProbesFromFilterOverlap(leafdevMapAndFilterTAIR9,chromosome=1,strand="forward",minPos=10000,maxPos=15000)
  sel <- tt$selection
  length(sel)
  head(sel)
  selfil <- tt$selectionFiltered
  length(selfil)
  head(selfil)
}

\keyword{manip}
