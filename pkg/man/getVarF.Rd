\name{getVarF}
\alias{getVarF}

\title{Get variance of model fit}

\description{Extract the variance of the fitted functional effects in the original data space.}

\usage{getVarF(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getVarF(leafdevFit)
  dim(tt)
  tt[,10:15]
}
}

\keyword{manip}
