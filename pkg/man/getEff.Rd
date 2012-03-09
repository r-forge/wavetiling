\name{getEff}
\alias{getEff}

\title{Get effects}

\description{Extract the estimated effects or contrasts of the wavelet-based transcriptome analysis.}

\usage{getEff(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getEff(leafdevInfCompare)
  dim(tt)
  tt[1:10,35:42]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
