\name{getNoProbes}
\alias{getNoProbes}

\title{Get number of probes}

\description{Extract the number of probes from a \code{WfmFit}-class object.}

\usage{getNoProbes(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{integer indicating the number of probes}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getNoProbes(leafdevFit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
