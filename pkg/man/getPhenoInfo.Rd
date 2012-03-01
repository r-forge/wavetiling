\name{getPhenoInfo}
\alias{getPhenoInfo}

\title{Get phenotypic information}

\description{Extract the phenotypic information for the tiling array experiment.}

\usage{getPhenoInfo(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{data frame}

\author{Kristof De Beuf <kristof.debeuf@ugent.be}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getPhenoInfo(leafdevFit)
}
}

\keyword{manip}
