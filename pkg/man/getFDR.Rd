\name{getFDR}
\alias{getFDR}

\title{get FDR}

\description{Extract the FDR for each test in the wavelet-based transcriptome analysis.}

\usage{getFDR(object)}

\arguments{
  \item{object}{object of class \code{WfmFit}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getFDR(leafdevInfCompare)
  dim(tt)
  tt[1:10,35:42]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
