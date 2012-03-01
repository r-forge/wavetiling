\name{getNoGroups}
\alias{getNoGroups}

\title{Extract the number of groups or samples}

\description{Accessor for the number of groups or samples in the tiling array experiment from a \code{WaveTilingFeatureSet} object}

\usage{getNoGroups(object)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
}

\value{numeric indicating the number of groups or samples in the tiling array experiment}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  leafdev <- addPheno(leafdev,noGroups=6,groupNames=c("day8","day9","day10","day11","day12","day13"),replics=rep(3,6))
  getNoGroups(leafdev)
}
}

\keyword{manip}
