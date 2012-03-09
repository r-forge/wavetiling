\name{getF}
\alias{getF}

\title{Get model fit}

\description{Extract the fitted functional effects in the original data space}

\usage{getF(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getF(leafdevFit)
  dim(tt)
  tt[,10:15]
}
}

\keyword{manip}
