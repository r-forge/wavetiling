\name{WfmInf-class}
\Rdversion{1.1}
\docType{class}
\alias{WfmInf-class}
\alias{WfmInf}
\alias{getAlpha}
\alias{getAlpha,WfmInf-method}
\alias{getDelta}
\alias{getDelta,WfmInf-method}
\alias{getEff}
\alias{getEff,WfmInf-method}
\alias{getFDR}
\alias{getFDR,WfmInf-method}
\alias{getGenomicRegions}
\alias{getGenomicRegions,WfmInf-method}
\alias{getRegions}
\alias{getRegions,WfmInf-method}
\alias{getSigProbes}
\alias{getSigProbes,WfmInf-method}
\alias{getTwoSided}
\alias{getTwoSided,WfmInf-method}
\alias{getSigGenes,WfmFit,WfmInf-method}
\alias{getNonAnnotatedRegions,WfmFit,WfmInf-method}
\alias{getVarEff}
\alias{getVarEff,WfmInf-method}
\alias{getGenomeInfo,WfmInf-method}
\alias{initialize,WfmInf-method}
\alias{plotWfm,WfmFit,WfmInf-method}
\alias{show,WfmInf-method}

\title{Class "WfmInf"}

\description{class to store outputs from the inference in the wavelet-based transcriptome analysis.}

\usage{

## Accessors

getAlpha(object)
getDelta(object)
getTwoSided(object)
getSigProbes(object)
getRegions(object)
getGenomicRegions(object)
getFDR(object)
getEff(object)
getVarEff(object)
}

\arguments{
  \item{object}{An instance of \code{WfmInf}-class.}
}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("WfmInf",  alpha, delta, two.sided, sigProbes, regions, GlocRegions, FDR, CI, eff, varEff)}.}

\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{delta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{two.sided}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sigProbes}:}{Object of class \code{"list"} ~~ }
    \item{\code{regions}:}{Object of class \code{"list"} ~~ }
    \item{\code{GlocRegions}:}{Object of class \code{"list"} ~~ }
    \item{\code{FDR}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{CI}:}{Object of class \code{"array"} ~~ }
    \item{\code{eff}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{varEff}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genome.info}:}{Object of class \code{"genomeInfo"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getAlpha}{\code{signature(object = "WfmInf")}: ... }
    \item{getDelta}{\code{signature(object = "WfmInf")}: ... }
    \item{getTwoSided}{\code{signature(object = "WfmInf")}: ... }
    \item{getSigProbes}{\code{signature(object = "WfmInf")}: ... }
    \item{getRegions}{\code{signature(object = "WfmInf")}: ... }
    \item{getGenomicRegions}{\code{signature(object = "WfmInf")}: ... }
    \item{getFDR}{\code{signature(object = "WfmInf")}: ... }
    \item{getEff}{\code{signature(object = "WfmInf")}: ... }
    \item{getVarEff}{\code{signature(object = "WfmInf")}: ... }
    \item{getGenomeInfo}{\code{signature(object = "WfmFit")}: ... }
    \item{initialize}{\code{signature(.Object = "WfmInf")}: ... }
    \item{show}{\code{signature(object = "WfmInf")}: ... }
    \item{plotWfm}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
    \item{getSigGenes}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
    \item{getNonAnnotatedRegions}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
	 }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a WfmInf object.
  \describe{
    \item{}{
      \code{getAlpha(x): Extract the alpha level of significance used in the wavelet-based analysis.}
    }
    \item{}{
      \code{getDelta(x): Extract the threshold values used in the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getTwoSided(x): Extract the direction of inference conducted in the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getSigProbes(x): Extract the significant probe ids for the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getRegions(x): Extract the significant regions from the wavelet-based transcriptome analysis. Regions are given in terms of the probe ids they map onto.}
    }
    \item{}{
      \code{getGenomicRegions(x): Extract the significant genomic regions from the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getFDR(x): Extract the FDR for each test in the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getEff(x): Extract the estimated effects or contrasts of the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getVarEff(x): Extract the estimated variances of the effects or contrasts in the wavelet-based transcriptome analysis.}
    }
    \item{}{
      \code{getGenomeInfo(x): Extract the genomic information.}
    }
  }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
  showClass("WfmInf")
  
  library(waveTilingData)
  data(leafdevInfCompare)
  tt1 <- getAlpha(leafdevInfCompare)
  tt2 <- getDelta(leafdevInfCompare)
  tt3 <- getTwoSided(leafdevInfCompare)
  tt4 <- getSigProbes(leafdevInfCompare)
  tt5 <- getRegions(leafdevInfCompare)
  tt6 <- getGenomicRegions(leafdevInfCompare)
  tt7 <- getFDR(leafdevInfCompare)
  tt8 <- getEff(leafdevInfCompare)
  tt9 <- getVarEff(leafdevInfCompare)
  tt10 <- getGenomeInfo(leafdevInfCompare)
}

\keyword{classes}


