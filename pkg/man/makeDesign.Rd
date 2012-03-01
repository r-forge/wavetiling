\name{makeDesign}
\alias{makeDesign}

\title{Construct design matrix}

\description{Helper function to construct a design matrix to be used in the wavelet-based transcriptome analysis.}

\usage{makeDesign(design=c("time","circadian","group","factorial"), replics, noGroups, factor.levels=NULL)}

\arguments{
  \item{design}{character indicating the design of the tiling array experiment. Currently, the following designs are implemented: \code{time} for a time-course design based on polynomial contrasts; \code{circadian} for circadian rhythm analysis; \code{group} for unordered one-factor designs; \code{factorial} for two-factor designs}
  \item{replics}{Numeric vector containing the number of replicates for each group. The vector length should be equal to the indicated number of groups.}
  \item{noGroups}{Number of groups in the tiling array experiment}
  \item{factor.levels}{Factor levels to use if applying two-factor design}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  makeDesign(design="time",replics=rep(3,6),noGroups=6)
}

\keyword{manip}
