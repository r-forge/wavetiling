\name{getDataOrigSpace}
\alias{getDataOrigSpace}

\title{Get data in original space}

\description{Extract the raw expression data in the original data space.}

\usage{getDataOrigSpace(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  tt <- getDataOrigSpace(leafdevFit)
  dim(tt)
  tt[1:10,10:15]
}
}

\keyword{manip}
