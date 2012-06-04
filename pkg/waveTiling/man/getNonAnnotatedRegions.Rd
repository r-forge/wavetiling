\name{getNonAnnotatedRegions}
\alias{getNonAnnotatedRegions}

\title{Get non-annotated regions}

\description{Extract the significant regions found in the wavelet-based transcriptome analysis that don't show any overlap with the existing annotation.}

\usage{getNonAnnotatedRegions(fit,inf,biomartObj)}

\arguments{
  \item{fit}{object of class \code{WfmFit}}
  \item{inf}{object of class \code{WfmInf}}
  \item{biomartObj}{object of class TranscriptDb representing an annotation database generated from BioMart.}
}

\value{GRangesList object with the non-annotated regions. The first element gives the regions with no annotation overlap on the strand used in the analysis, the second element gives the regions with no annotation overlap on both strands.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  library(waveTilingData)
  library(TxDb.Athaliana.BioMart.plantsmart12)
  data(leafdevFit)
  data(leafdevInfCompare)
  nonAnnoCompare <- getNonAnnotatedRegions(fit=leafdevFit,inf=leafdevInfCompare,	biomartObj=TxDb.Athaliana.BioMart.plantsmart12)
  nonAnnoCompare
}

\keyword{manip}
