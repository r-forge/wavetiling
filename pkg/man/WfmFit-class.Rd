\name{WfmFit-class}
\Rdversion{1.1}
\docType{class}
\alias{WfmFit-class}
\alias{show,WfmFit-method}
\alias{initialize,WfmFit-method}
\alias{getProbePosition,WfmFit-method}
\alias{getNoProbes,WfmFit-method}
\alias{getBetaWav,WfmFit-method}
\alias{getVarBetaWav,WfmFit-method}
\alias{getSmoothPar,WfmFit-method}
\alias{getVarEps,WfmFit-method}
\alias{getGenomeInfo,WfmFit-method}
\alias{getChromosome,WfmFit-method}
\alias{getStrand,WfmFit-method}
\alias{getMinPos,WfmFit-method}
\alias{getMaxPos,WfmFit-method}
\alias{getNoLevels,WfmFit-method}
\alias{getDesignMatrix,WfmFit-method}
\alias{getPhenoInfo,WfmFit-method}
\alias{getDataOrigSpace,WfmFit-method}
\alias{getDataWaveletSpace,WfmFit-method}
\alias{getWaveletFilter,WfmFit-method}
\alias{getKj,WfmFit-method}
\alias{getPrior,WfmFit-method}
\alias{getF,WfmFit-method}
\alias{getVarF,WfmFit-method}
\alias{wfm.inference,WfmFit-method}

\title{Class "WfmFit"}
\description{class to store model fits in the wavelet-based transcriptome analysis.}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("WfmFit", betaWav, varbetaWav, smoothPar, varEps, dataOrigSpace, dataWaveletSpace, design.matrix, phenoData, genome.info, n.levels, probePosition, wave.filt, Kj, prior, F, varF, P, Z, noGroups, replics)}.}

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

\section{Methods}{
  \describe{
    \item{getBetaWav}{\code{signature(object = "WfmFit")}: ... }
    \item{getChromosome}{\code{signature(object = "WfmFit")}: ... }
    \item{getDataOrigSpace}{\code{signature(object = "WfmFit")}: ... }
    \item{getDataWaveletSpace}{\code{signature(object = "WfmFit")}: ... }
    \item{getDesignMatrix}{\code{signature(object = "WfmFit")}: ... }
    \item{getF}{\code{signature(object = "WfmFit")}: ... }
    \item{getGenomeInfo}{\code{signature(object = "WfmFit")}: ... }
    \item{getKj}{\code{signature(object = "WfmFit")}: ... }
    \item{getMaxPos}{\code{signature(object = "WfmFit")}: ... }
    \item{getMinPos}{\code{signature(object = "WfmFit")}: ... }
    \item{getNoLevels}{\code{signature(object = "WfmFit")}: ... }
    \item{getNoProbes}{\code{signature(object = "WfmFit")}: ... }
    \item{getPhenoInfo}{\code{signature(object = "WfmFit")}: ... }
    \item{getPrior}{\code{signature(object = "WfmFit")}: ... }
    \item{getProbePosition}{\code{signature(object = "WfmFit")}: ... }
    \item{getSmoothPar}{\code{signature(object = "WfmFit")}: ... }
    \item{getStrand}{\code{signature(object = "WfmFit")}: ... }
    \item{getVarBetaWav}{\code{signature(object = "WfmFit")}: ... }
    \item{getVarEps}{\code{signature(object = "WfmFit")}: ... }
    \item{getVarF}{\code{signature(object = "WfmFit")}: ... }
    \item{getWaveletFilter}{\code{signature(object = "WfmFit")}: ... }
    \item{initialize}{\code{signature(.Object = "WfmFit")}: ... }
    \item{show}{\code{signature(object = "WfmFit")}: ... }
    \item{wfm.inference}{\code{signature(object = "WfmFit")}: ... }
    }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("WfmFit")
}

\keyword{classes}
