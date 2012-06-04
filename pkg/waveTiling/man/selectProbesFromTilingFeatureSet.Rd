\name{selectProbesFromTilingFeatureSet}
\alias{selectProbesFromTilingFeatureSet}

\title{select probes from \code{WaveTilingFeatureSet} object}

\description{Extract the probe indices from a \code{WaveTilingFeatureSet} object that map to a region between two specified genomic positions}

\usage{selectProbesFromTilingFeatureSet(object, chromosome, strand=c("forward","reverse"), minPos, maxPos)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
  \item{chromosome}{chromosome}
  \item{strand}{strand}
  \item{minPos}{minimum genomic position}
  \item{maxPos}{maximum genomic position}
}

\value{vector of integers indicating the probe indices}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  \dontrun{
    data(leafdevBQ)
    tt <- selectProbesFromTilingFeatureSet(leafdevBQ,chromosome=1,strand="forward",minPos=22000,maxPos=26000)
    length(tt)
    head(tt)
  }
}

\keyword{manip}
