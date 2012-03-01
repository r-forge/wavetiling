\name{getMinPos}
\alias{getMinPos}

\title{Get minimum genomic position}

\description{Extract the minimum genomic probe position from a \code{WfmFit}-class object}

\usage{getMinPos(object)}

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
