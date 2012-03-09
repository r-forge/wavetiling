\name{getDataWaveletSpace}
\alias{getDataWaveletSpace}

\title{Get data in the wavelet space}

\description{Extract the raw data in the wavelet space, i.e. the wavelet coefficients.
}

\usage{getDataWaveletSpace(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getDataWaveletSpace(leafdevFit)
  dim(tt)
  tt[1:10,10:15]
}
}

\keyword{manip}
