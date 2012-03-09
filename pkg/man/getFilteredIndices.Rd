\name{getFilteredIndices}
\alias{getFilteredIndices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{get filtered indices}

\description{Extract the filtered probe indices from a \code{MapFilterProbe}-object.}

\usage{getFilteredIndices(object)}

\arguments{
  \item{object}{object of class \code{MapFilterProbe}}
}

\value{numeric vector containing the filtered probe indices}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevMapAndFilterTAIR9)
  tt <- getFilteredIndices(leafdevMapAndFilterTAIR9)
  length(tt)
  head(tt)
}
}

\keyword{manip}
