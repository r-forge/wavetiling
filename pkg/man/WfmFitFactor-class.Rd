\name{WfmFitFactor-class}
\Rdversion{1.1}
\docType{class}
\alias{WfmFitFactor-class}
\alias{show,WfmFitFactor-method}
\alias{initialize,WfmFitFactor-method}

\title{Class "WfmFitFactor"}
\description{class to store model fits with factorial design in the wavelet-based transcriptome analysis.}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("WfmFitFactor")}.}

\section{Slots}{
  \describe{
    \item{\code{betaWav}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varbetaWav}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{smoothPar}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varEps}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dataOrigSpace}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{dataWaveletSpace}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{design.matrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{genome.info}:}{Object of class \code{"genomeInfo"} ~~ }
    \item{\code{n.levels}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{probePosition}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{wave.filt}:}{Object of class \code{"character"} ~~ }
    \item{\code{Kj}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{prior}:}{Object of class \code{"character"} ~~ }
    \item{\code{F}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varF}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{P}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Z}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{noGroups}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{replics}:}{Object of class \code{"numeric"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{WfmFit}"}, directly.
}

\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "WfmFitFactor")}: ... }
    \item{show}{\code{signature(object = "WfmFitFactor")}: ... }
    }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("WfmFitFactor")
}

\keyword{classes}
