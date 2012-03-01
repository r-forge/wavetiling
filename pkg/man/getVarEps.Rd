\name{getVarEps}
\alias{getVarEps}

\title{Get residual variance in wavelet space}

\description{Extract the estimated residual variance in the wavelet space. One variance parameter is estimated per wavelet level.}

\usage{getVarEps(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric vector}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getVarEps(leafdevFit)
}
}

\keyword{manip}
