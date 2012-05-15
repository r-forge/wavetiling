\name{MapFilterProbe-class}
\Rdversion{1.1}
\docType{class}
\alias{MapFilterProbe-class}
\alias{MapFilterProbe}
\alias{getChromosome,MapFilterProbe-method}
\alias{getChromosome}
\alias{getFilteredIndices,MapFilterProbe-method}
\alias{getFilteredIndices}
\alias{getPosition,MapFilterProbe-method}
\alias{getPosition}
\alias{getStrand,MapFilterProbe-method}
\alias{getStrand}
\alias{initialize,MapFilterProbe-method}
\alias{selectProbesFromFilterOverlap,MapFilterProbe-method}
\alias{show,MapFilterProbe-method}

\title{Class "MapFilterProbe"}

\description{class to store probe information after remapping and/or filtering of probes.}

\usage{

## Accessors

getChromosome(object)
getFilteredIndices(object)
getPosition(object)
getStrand(object)
}

\arguments{
  \item{object}{An instance of \code{MapFilterProbe}-class.}
}

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

\section{Accessors}{
  In the following code snippets, \code{x} is a MapFilterProbe object.
  \describe{
    \item{}{
      \code{getChromosome(x): Extract the chromosome identifiers.}
    }
    \item{}{
      \code{getFilteredIndices(x): Extract the filtered probe indices.}
    }
    \item{}{
      \code{getPosition(x): Extract the genomic position of the filtered probes.}
    }
    \item{}{
      \code{getStrand(x): Extract the strand orientation info for the filtered probes.}
    }
  }
}


\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  showClass("MapFilterProbe")
  \dontrun{
    library(waveTilingData)
    data(leafdevMapAndFilterTAIR9)
    tt1 <- getChromsome(leafdevMapAndFilterTAIR9)
    tt2 <- getFilteredIndices(leafdevMapAndFilterTAIR9)
    tt3 <- getPosition(leafdevMapAndFilterTAIR9)
    tt4 <- getStrand(leafdevMapAndFilterTAIR9)
  }
}

\keyword{classes}


