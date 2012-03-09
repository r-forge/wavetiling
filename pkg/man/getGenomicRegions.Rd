\name{getGenomicRegions}
\alias{getGenomicRegions}

\title{Get significantg genomic regions}

\description{Extract the significant genomic regions from the wavelet-based transcriptome analysis.}

\usage{getGenomicRegions(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{list of IRanges containing genomic start and end positions of the significant regions for the corresponding analysis.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getGenomicRegions(leafdevInfCompare)
  class(tt)
  class(tt[[2]])
  tt[[2]]
}
}

\keyword{manip}
