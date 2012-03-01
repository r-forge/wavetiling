\name{getGenomeInfo}
\alias{getGenomeInfo}

\title{Get genomic information}

\description{Extract the genomic information from a \code{WfmFit}-class object.}

\usage{getGenomeInfo(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{object of class genomeInfo which contains the chromosome number, the strand orientation and the minimum and the maximum genomic probe position.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getGenomeInfo(leafdevFit)
}
}

\keyword{manip}
