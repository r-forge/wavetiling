\name{getReplics}
\alias{getReplics}

\title{Get the number of replicates}

\description{Accessor for the number of replicates in the tiling array experiment from a \code{WaveTilingFeatureSet} object}

\usage{getReplics(object)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
}

\value{numeric vector with number of replicates for each of the groups or samples in the tiling array experiment}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  leafdev <- addPheno(leafdev,noGroups=6,groupNames=c("day8","day9","day10","day11","day12","day13"),replics=rep(3,6))
  getReplics(leafdev)
}
}

\keyword{manip}
