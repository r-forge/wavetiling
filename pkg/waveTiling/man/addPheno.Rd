\name{addPheno}
\alias{addPheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add phenotypic info to \code{WaveTilingFeatureSet}}

\description{Function to add phenotypic information such as sample names or group names to a \code{WaveTilingFeatureSet}-class object}

\usage{addPheno(object, noGroups, groupNames, replics, ...)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
  \item{noGroups}{Number of groups in the tiling array experiment}
  \item{groupNames}{Vector containing the group or sample names in the tiling array experiment. The vector length should be equal to the indicated number of groups.}
  \item{replics}{Numeric vector containing the number of replicates for each group. The vector length should be equal to the indicated number of groups.}
  \item{\dots}{other arguments}
}

\value{object of class \code{WaveTilingFeatureSet} annotated with the phenotypic data}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  leafdev <- addPheno(leafdev,noGroups=6,groupNames=c("day8","day9","day10","day11","day12","day13"),replics=rep(3,6))
  leafdev
}
}

\keyword{manip}

