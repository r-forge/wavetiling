\name{getMaxPos}
\alias{getMaxPos}

\title{Get maximum genomic position}

\description{Extract the minimum genomic probe position from a \code{WfmFit}-class object}

\usage{getMaxPos(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{integer}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getMinPos(leafdevFit)
}
}

\keyword{manip}
