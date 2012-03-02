\name{makeContrasts}
\alias{makeContrasts}

\title{Construct contrast matrix}

\description{Helper function to construct a contrast matrix to be used in the inference procedure of the wavelet-based transcriptome analysis when conducting a pairwise comparison analysis.}

\usage{makeContrasts(contrasts, nlevels)}

\arguments{
  \item{contrasts}{\code{compare}: contrasts for pairwise comparison analysis.}
  \item{nlevels}{Number of groups for pairwise comparison analysis.}
}

\value{numeric matrix}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  makeContrasts(contrasts="compare",nlevels=6)
}

\keyword{manip}
