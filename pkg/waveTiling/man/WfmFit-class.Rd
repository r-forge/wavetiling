\name{WfmFit-class}
\Rdversion{1.1}
\docType{class}
\alias{WfmFit-class}
\alias{WfmFit}
\alias{show,WfmFit-method}
\alias{initialize,WfmFit-method}
\alias{getProbePosition}
\alias{getProbePosition,WfmFit-method}
\alias{getNoProbes}
\alias{getNoProbes,WfmFit-method}
\alias{getBetaWav}
\alias{getBetaWav,WfmFit-method}
\alias{getVarBetaWav}
\alias{getVarBetaWav,WfmFit-method}
\alias{getSmoothPar}
\alias{getSmoothPar,WfmFit-method}
\alias{getVarEps}
\alias{getVarEps,WfmFit-method}
\alias{getGenomeInfo}
\alias{getGenomeInfo,WfmFit-method}
\alias{getChromosome,WfmFit-method}
\alias{getStrand,WfmFit-method}
\alias{getMinPos}
\alias{getMinPos,WfmFit-method}
\alias{getMaxPos}
\alias{getMaxPos,WfmFit-method}
\alias{getNoLevels}
\alias{getNoLevels,WfmFit-method}
\alias{getDesignMatrix}
\alias{getDesignMatrix,WfmFit-method}
\alias{getPhenoInfo}
\alias{getPhenoInfo,WfmFit-method}
\alias{getDataOrigSpace}
\alias{getDataOrigSpace,WfmFit-method}
\alias{getDataWaveletSpace}
\alias{getDataWaveletSpace,WfmFit-method}
\alias{getWaveletFilter}
\alias{getWaveletFilter,WfmFit-method}
\alias{getKj}
\alias{getKj,WfmFit-method}
\alias{getPrior}
\alias{getPrior,WfmFit-method}
\alias{getF}
\alias{getF,WfmFit-method}
\alias{getVarF}
\alias{getVarF,WfmFit-method}
\alias{wfm.inference,WfmFit-method}

\title{Class "WfmFit"}

\description{class to store model fits in the wavelet-based transcriptome analysis.}

\usage{

## Accessors

getProbePosition(object)
getNoProbes(object)
getBetaWav(object)
getVarBetaWav(object)
getSmoothPar(object)
getVarEps(object)
getGenomeInfo(object)
getMinPos(object)
getMaxPos(object)
getNoLevels(object)
getDesignMatrix(object)
getPhenoInfo(object)
getDataOrigSpace(object)
getDataWaveletSpace(object)
getWaveletFilter(object)
getKj(object)
getPrior(object)
getF(object)
getVarF(object)
}

\arguments{
  \item{object}{An instance of \code{WfmFit}-class.}
}

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

\section{Accessors}{
  In the following code snippets, \code{x} is a WfmFit object.
  \describe{
    \item{}{
      \code{getBetaWav(x): Extract the fitted effect functions in the wavelet space.}
    }
    \item{}{
      \code{getChromsome(x): Extract the chromosome identifiers.}
    }
    \item{}{
      \code{getDataOrigSpace(x): Extract the raw expression data in the original data space.}
    }
    \item{}{
      \code{getDataWaveletSpace(x): Extract the raw data in the wavelet space, i.e. the wavelet coefficients.}
    }
    \item{}{
      \code{getDesignMatrix(x): Extract the design matrix used in the wavelet-based analysis.}
    }
    \item{}{
      \code{getF(x): Extract the fitted functional effects in the original data space.}
    }
    \item{}{
      \code{getGenomeInfo(x): Extract the genomic information.}
    }
    \item{}{
      \code{getKj(x): Extract the number of wavelet coefficients estimated per wavelet level.}
    }
    \item{}{
      \code{getMaxPos(x): Extract the maximum genomic probe position.}
    }
    \item{}{
      \code{getMinPos(x): Extract the minimum genomic probe position.}
    }
    \item{}{
      \code{getNoLevels(x): Extract the number of levels in in the wavelet decomposition when fitting the wavelet-based functional model.}
    }
    \item{}{
      \code{getNoProbes(x): Extract the number of probes.}
    }
    \item{}{
      \code{getPhenoInfo(x): Extract the phenotypic information for the tiling array experiment.}
    }
    \item{}{
      \code{getPrior(x): Extract the the type or distribution of the prior imposed on the functional effects in the wavelet space.}
    }
    \item{}{
      \code{getProbePosition(x): Extract probe position.}
    }
    \item{}{
      \code{getSmoothPar(x): Extract the estimated smoothing parameters that control the regularization of the effect functions in the wavelet space.}
    }
    \item{}{
      \code{getStrand(x): Extract the strand orientation info.}
    }
    \item{}{
      \code{getVarBetaWav(x): Extract the variance of the fitted effect functions in the wavelet space.}
    }
    \item{}{
      \code{getVarEps(x): Extract the estimated residual variance in the wavelet space. One variance parameter is estimated per wavelet level.}
    }
    \item{}{
      \code{getVarF(x): Extract the variance of the fitted functional effects in the original data space.}
    }
    \item{}{
      \code{getWaveletFilter(x): Extract the wavelet filter used to transform the data from the original space to the wavelet space.}
    }
  }
}


\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  showClass("WfmFit")
  
  library(waveTilingData)
  data(leafdevFit)
  tt1 <- getBetaWav(leafdevFit)
  tt2 <- getChromosome(leafdevFit)
  tt3 <- getDataOrigSpace(leafdevFit)
  tt4 <- getDataWaveletSpace(leafdevFit)
  tt5 <- getDesignMatrix(leafdevFit)
  tt6 <- getF(leafdevFit)
  tt7 <- getGenomeInfo(leafdevFit)
  tt8 <- getKj(leafdevFit)
  tt9 <- getMaxPos(leafdevFit)
  tt10 <- getMinPos(leafdevFit)
  tt11 <- getNoLevels(leafdevFit)
  tt12 <- getNoProbes(leafdevFit)
  tt13 <- getPhenoInfo(leafdevFit)
  tt14 <- getPrior(leafdevFit)
  tt15 <- getProbePosition(leafdevFit)
  tt16 <- getSmoothPar(leafdevFit)
  tt17 <- getStrand(leafdevFit)
  tt18 <- getVarBetaWav(leafdevFit)
  tt19 <- getVarEps(leafdevFit)
  tt20 <- getVarF(leafdevFit)
  tt21 <- getWaveletFilter(leafdevFit)
}

\keyword{classes}
