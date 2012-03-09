\name{getKj}
\alias{getKj}

\title{Get Kj}

\description{Extract the number of wavelet coefficients estimated per wavelet level.}

\usage{getKj(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric vector}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getKj(leafdevFit)
}
}

\keyword{manip}
