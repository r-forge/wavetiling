\name{getVarEff}
\alias{getVarEff}

\title{get variance of effects}

\description{Extract the estimated variances of the effects or contrasts in the wavelet-based transcriptome analysis.}

\usage{getVarEff(object)}

\arguments{
  \item{object}{object of class \code{WfmInf}}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
\dontrun{
  library(waveTilingData)
  data(leafdevInfCompare)
  tt <- getVarEff(leafdevInfCompare)
  dim(tt)
  tt[1:10,35:42]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
