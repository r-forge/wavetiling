\name{getProbePosition}
\alias{getProbePosition}

\title{get probe position}

\description{Extract probe position from \code{WfmFit}-class object.}

\usage{getProbePosition(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{vector of integers containing the genomic probe positions}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}


\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getProbePosition(leafdevFit)
  length(tt)
  head(tt)
}
}

\keyword{manip}
