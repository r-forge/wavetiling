\name{getChromosome}
\alias{getChromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get chromosome identifiers}

\description{Extract the chromosome identifiers from a \code{MapFilterProbe}- or a \code{WfmFit}-object.}

\usage{getChromosome(object)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{MapFilterProbe} or of class \code{WfmFit}}
}

\value{If applied on a \code{MapFilterProbe} object, a vector of characters is returned, providing the chromosome identification of each filtered element. If applied on a \code{WfmFit} object, a numeric is returned providing the chromsome number for which the model for transcriptome analysis was fitted.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevMapAndFilterTAIR9)
  tt <- getChromosome(leafdevMapAndFilterTAIR9)
  head(tt)
  
  data(leafdevFit)
  getChromosome(leafdevFit)
}
}

\keyword{manip}
