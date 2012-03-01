\name{getStrand}
\alias{getStrand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get chromosome strand information}

\description{Extract the strand orientation info for the filtered probes of a \code{MapFilterProbe}-object or from a \code{WfmFit}-object.}

\usage{getStrand(object)}

\arguments{
  \item{object}{object of class \code{MapFilterProbe} or of class \code{WfmFit}}
}

\value{If applied on a \code{MapFilterProbe} object, a vector of characters is returned, providing the chromosome strand orientation of each filtered element. If applied on a \code{WfmFit} object, a character is returned providing the chromsome strand orientation for which the model for transcriptome analysis was fitted.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevMapAndFilterTAIR9)
  tt <- getStrand(leafdevMapAndFilterTAIR9)
  head(tt)
  
  data(leafdevFit)
  getStrand(leafdevFit)
}
}

\keyword{manip}
