\name{getWaveletFilter}
\alias{getWaveletFilter}

\title{Get wavelet filter used}

\description{Extract the wavelet filter used to transform the data from the original space to the wavelet space.}

\usage{getWaveletFilter(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{character}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevFit)
  getWaveletFilter(leafdevFit)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

