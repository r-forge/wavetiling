\name{Wfm-class}
\Rdversion{1.1}
\docType{class}
\alias{Wfm-class}
\alias{getAlpha,Wfm-method}
\alias{getBetaMAP,Wfm-method}
\alias{getChromosome,Wfm-method}
\alias{getCI,Wfm-method}
\alias{getDataOrigSpace,Wfm-method}
\alias{getDataWaveletSpace,Wfm-method}
\alias{getDelta,Wfm-method}
\alias{getDesign,Wfm-method}
\alias{getEff,Wfm-method}
\alias{getFDR,Wfm-method}
\alias{getF,Wfm-method}
\alias{getGenomeInfo,Wfm-method}
\alias{getGenomicRegions,Wfm-method}
\alias{getKj,Wfm-method}
\alias{getMaxPos,Wfm-method}
\alias{getMinPos,Wfm-method}
\alias{getNoLevels,Wfm-method}
\alias{getNoProbes,Wfm-method}
\alias{getPhenoInfo,Wfm-method}
\alias{getPrior,Wfm-method}
\alias{getProbePosition,Wfm-method}
\alias{getRegions,Wfm-method}
\alias{getRescale,Wfm-method}
\alias{getSigProbes,Wfm-method}
\alias{getSmoothPar,Wfm-method}
\alias{getStrand,Wfm-method}
\alias{getTwoSided,Wfm-method}
\alias{getVarBetaMAP,Wfm-method}
\alias{getVarEff,Wfm-method}
\alias{getVarEps,Wfm-method}
\alias{getVarF,Wfm-method}
\alias{getWaveletFilter,Wfm-method}
\alias{getWfmMethod,Wfm-method}
\alias{initialize,Wfm-method}
\alias{plotWfm,Wfm-method}
\alias{show,Wfm-method}

\title{Class "Wfm"}
\description{
##%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Wfm", betaMAP, varbetaMAP, smoothPar, varEps, dataOrigSpace, dataWaveletSpace, design, phenoData, method, genome.info, n.levels, probePosition, wave.filt, Kj, prior, alpha, delta, two.sided, rescale, sigProbes, regions, GlocRegions, FDR, CI, F, varF, eff, varEff)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{betaMAP}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varbetaMAP}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{smoothPar}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varEps}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dataOrigSpace}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{dataWaveletSpace}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{design}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{genome.info}:}{Object of class \code{"genomeInfo"} ~~ }
    \item{\code{n.levels}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{probePosition}:}{Object of class \code{"vector"} ~~ }
    \item{\code{wave.filt}:}{Object of class \code{"character"} ~~ }
    \item{\code{Kj}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{prior}:}{Object of class \code{"character"} ~~ }
    \item{\code{alpha}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{delta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{two.sided}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rescale}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{sigProbes}:}{Object of class \code{"list"} ~~ }
    \item{\code{regions}:}{Object of class \code{"list"} ~~ }
    \item{\code{GlocRegions}:}{Object of class \code{"list"} ~~ }
    \item{\code{FDR}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{CI}:}{Object of class \code{"array"} ~~ }
    \item{\code{F}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varF}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{eff}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varEff}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getAlpha}{\code{signature(object = "Wfm")}: ... }
    \item{getBetaMAP}{\code{signature(object = "Wfm")}: ... }
    \item{getChromosome}{\code{signature(object = "Wfm")}: ... }
    \item{getCI}{\code{signature(object = "Wfm")}: ... }
    \item{getDataOrigSpace}{\code{signature(object = "Wfm")}: ... }
    \item{getDataWaveletSpace}{\code{signature(object = "Wfm")}: ... }
    \item{getDelta}{\code{signature(object = "Wfm")}: ... }
    \item{getDesign}{\code{signature(object = "Wfm")}: ... }
    \item{getEff}{\code{signature(object = "Wfm")}: ... }
    \item{getFDR}{\code{signature(object = "Wfm")}: ... }
    \item{getF}{\code{signature(object = "Wfm")}: ... }
    \item{getGenomeInfo}{\code{signature(object = "Wfm")}: ... }
    \item{getGenomicRegions}{\code{signature(object = "Wfm")}: ... }
    \item{getKj}{\code{signature(object = "Wfm")}: ... }
    \item{getMaxPos}{\code{signature(object = "Wfm")}: ... }
    \item{getMinPos}{\code{signature(object = "Wfm")}: ... }
    \item{getNoLevels}{\code{signature(object = "Wfm")}: ... }
    \item{getNoProbes}{\code{signature(object = "Wfm")}: ... }
    \item{getPhenoInfo}{\code{signature(object = "Wfm")}: ... }
    \item{getPrior}{\code{signature(object = "Wfm")}: ... }
    \item{getProbePosition}{\code{signature(object = "Wfm")}: ... }
    \item{getRegions}{\code{signature(object = "Wfm")}: ... }
    \item{getRescale}{\code{signature(object = "Wfm")}: ... }
    \item{getSigProbes}{\code{signature(object = "Wfm")}: ... }
    \item{getSmoothPar}{\code{signature(object = "Wfm")}: ... }
    \item{getStrand}{\code{signature(object = "Wfm")}: ... }
    \item{getTwoSided}{\code{signature(object = "Wfm")}: ... }
    \item{getVarBetaMAP}{\code{signature(object = "Wfm")}: ... }
    \item{getVarEff}{\code{signature(object = "Wfm")}: ... }
    \item{getVarEps}{\code{signature(object = "Wfm")}: ... }
    \item{getVarF}{\code{signature(object = "Wfm")}: ... }
    \item{getWaveletFilter}{\code{signature(object = "Wfm")}: ... }
    \item{getWfmMethod}{\code{signature(object = "Wfm")}: ... }
    \item{initialize}{\code{signature(.Object = "Wfm")}: ... }
    \item{plotWfm}{\code{signature(object = "Wfm")}: ... }
    \item{show}{\code{signature(object = "Wfm")}: ... }
	 }
}
\references{
## %%  ~~put references to the literature/web site here~~
}
\author{
Kristof De Beuf <kristof.debeuf@ugent.be>
}
\note{
## %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
## %%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Wfm")
}
\keyword{classes}
