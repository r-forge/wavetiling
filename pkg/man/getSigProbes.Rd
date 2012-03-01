\name{getSigProbes}
\alias{getSigProbes}

\title{Get significant probes}

\description{Extract the significant probe ids for the wavelet-based transcriptome analysis.}

\usage{getSigProbes(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{list of numeric vectors. Each element in the list contains a vector with the probe ids of the significant probes for the corresponding analysis.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getSigProbes(leafdevInfCompare)
  class(tt)
  head(tt[[2]])
  length(tt[[2]])
}
}

\keyword{manip}
