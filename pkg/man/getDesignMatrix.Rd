\name{getDesignMatrix}
\alias{getDesignMatrix}

\title{Get design matrix}

\description{Extract the design matrix used in the wavelet-based analysis.}

\usage{getDesignMatrix(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getDesignMatrix(leafdevFit)
}
}

\keyword{manip}
