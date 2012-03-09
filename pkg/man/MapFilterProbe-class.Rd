\name{MapFilterProbe-class}
\Rdversion{1.1}
\docType{class}
\alias{MapFilterProbe-class}
\alias{getChromosome,MapFilterProbe-method}
\alias{getFilteredIndices,MapFilterProbe-method}
\alias{getPosition,MapFilterProbe-method}
\alias{getStrand,MapFilterProbe-method}
\alias{initialize,MapFilterProbe-method}
\alias{selectProbesFromFilterOverlap,MapFilterProbe-method}
\alias{show,MapFilterProbe-method}

\title{Class "MapFilterProbe"}

\description{class to store probe information after remapping and/or filtering of probes.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mapFilterProbe", filteredIndices, chromosome, position, strand)}.}

\section{Slots}{
  \describe{
    \item{\code{filteredIndices}:}{Object of class \code{"vector"} ~~ }
    \item{\code{chromosome}:}{Object of class \code{"vector"} ~~ }
    \item{\code{position}:}{Object of class \code{"vector"} ~~ }
    \item{\code{strand}:}{Object of class \code{"vector"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getChromosome}{\code{signature(object = "MapFilterProbe")}: ... }
    \item{getFilteredIndices}{\code{signature(object = "MapFilterProbe")}: ... }
    \item{getPosition}{\code{signature(object = "MapFilterProbe")}: ... }
    \item{getStrand}{\code{signature(object = "MapFilterProbe")}: ... }
    \item{initialize}{\code{signature(.Object = "MapFilterProbe")}: ... }
    \item{selectProbesFromFilterOverlap}{\code{signature(object = "MapFilterProbe")}: ... }
    \item{show}{\code{signature(object = "MapFilterProbe")}: ... }
	 }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("MapFilterProbe")
}

\keyword{classes}


