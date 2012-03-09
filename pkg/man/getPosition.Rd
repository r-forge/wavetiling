\name{getPosition}
\alias{getPosition}

\title{get probe positions}

\description{Extract the genomic position of the filtered probes from a\code{MapFilterProbe}-object.}

\usage{getPosition(object)}

\arguments{
  \item{object}{object of class \code{MapFilterProbe}}
}

\value{numeric vector containing the genomic positions of the filtered probes}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevMapAndFilterTAIR9)
  tt <- getPosition(leafdevMapAndFilterTAIR9)
  length(tt)
  head(tt)
}
}

\keyword{manip}
