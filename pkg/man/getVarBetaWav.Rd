\name{getVarBetaWav}
\alias{getVarBetaWav}

\title{Get variance of wavelet space effect functions}

\description{Extract the variance of the fitted effect functions in the wavelet space}

\usage{getVarBetaWav(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getVarBetaWav(leafdevFit)
  dim(tt)
  tt[,10:15]
}
}

\keyword{manip}
