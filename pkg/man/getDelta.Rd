\name{getDelta}
\alias{getDelta}

\title{Get threshold values delta}

\description{Extract the threshold values used in the wavelet-based transcriptome analysis.}

\usage{getDelta(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric vector with the threshold values used for each type of test}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  getDelta(leafdevInfCompare)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
