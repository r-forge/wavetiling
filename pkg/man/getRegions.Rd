\name{getRegions}
\alias{getRegions}

\title{Get singificant regions}

\description{Extract the significant regions from the wavelet-based transcriptome analysis. Regions are given in terms of the probe ids they map onto.}

\usage{getRegions(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{list of IRanges containing start and end positions of the significant regions for the corresponding analysis.}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getRegions(leafdevInfCompare)
  class(tt)
  class(tt[[2]])
  tt[[2]]
}
}

\keyword{manip}
