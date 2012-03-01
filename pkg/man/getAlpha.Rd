\name{getAlpha}
\alias{getAlpha}

\title{Get significance level alpha}

\description{Accessor function to call the alpha level of significance used in the wavelet-based analysis.}

\usage{getAlpha(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{numeric vector with significance levels used for each type of test}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  getAlpha(leafdevInfCompare)
}
}

\keyword{manip}
