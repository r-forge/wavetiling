\name{getBetaWav}
\alias{getBetaWav}

\title{Get wavelet space effect functions}

\description{Extract the fitted effect functions in the wavelet space}

\usage{getBetaWav(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getBetaWav(leafdevFit)
  dim(tt)
  tt[,10:15]
}
}

\keyword{manip}
