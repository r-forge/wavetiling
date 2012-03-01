\name{getTwoSided}
\alias{getTwoSided}

\title{One- or two-sided inference}

\description{Extract the direction of inference conducted in the wavelet-based transcriptome analysis.}

\usage{getTwoSided(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{numeric vector indicating for each type of test whether a two-sided analysis was used (1) or not (0)}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  getTwoSided(leafdevInfCompare)
}
}

\keyword{manip}
