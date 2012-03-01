\name{getPrior}
\alias{getPrior}

\title{get prior type or distribution}

\description{Extract the the type or distribution of the prior imposed on the functional effects in the wavelet space.}

\usage{getPrior(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{character}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getPrior(leafdevFit)
}
}

\keyword{manip}
