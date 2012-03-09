\name{getNoLevels}
\alias{getNoLevels}

\title{Get number of levels in wavelet decomposition}

\description{Extract the number of levels in in the wavelet decomposition when fitting the wavelet-based functional model.}

\usage{getNoLevels(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{integer}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getNoLevels(leafdevFit)
}
}

\keyword{manip}
