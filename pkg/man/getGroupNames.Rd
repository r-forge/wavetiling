\name{getGroupNames}
\alias{getGroupNames}

\title{Get the group or sample names}
\description{Accessor for the group or sample names in the tiling array experiment from a \code{WaveTilingFeatureSet} object}

\usage{getGroupNames(object)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
}

\value{character vector containing the group or sample names in the tiling array experiment}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  leafdev <- addPheno(leafdev,noGroups=6,groupNames=c("day8","day9","day10","day11","day12","day13"),replics=rep(3,6))
  getGroupNames(leafdev)
}
}

\keyword{manip}
