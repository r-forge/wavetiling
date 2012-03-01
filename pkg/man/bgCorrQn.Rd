\name{bgCorrQn}
\alias{bgCorrQn}

\title{Background correction and quantile normalization}

\description{Function to perform background correction and quantile normalization of tiling arrays}

\usage{bgCorrQn(object, useMapFilter=NULL)}

\arguments{
  \item{object}{object of class \code{WaveTilingFeatureSet}}
  \item{useMapFilter}{NULL or object of class \code{mapFilterProbe} indicating the probes to use for background correction and quantile normalization}
  }

\value{object of class \code{WaveTilingFeatureSet} containing the background-corrected and quantile-normalized intensity signals}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  data(leafdev)
  leafdev <- as(leafdev,"WaveTilingFeatureSet")
  data(leafdevMapAndFilterTAIR9)
  leafdevBQ <- bgCorrQn(leafdev,useMapFilter=leafdevMapAndFilterTAIR9)
}
}

\keyword{manip}
