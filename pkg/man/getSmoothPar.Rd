\name{getSmoothPar}
\alias{getSmoothPar}

\title{Get smoothing parameters}

\description{Extract the estimated smoothing parameters that control the regularization of the effect functions in the wavelet space.}

\usage{getSmoothPar(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getSmoothPar(leafdevFit)
  dim(tt)
  tt[,10:15]
}
}

\keyword{manip}
