\name{WfmInf-class}
\Rdversion{1.1}
\docType{class}
\alias{WfmInf-class}
\alias{getAlpha,WfmInf-method}
\alias{getDelta,WfmInf-method}
\alias{getEff,WfmInf-method}
\alias{getFDR,WfmInf-method}
\alias{getGenomicRegions,WfmInf-method}
\alias{getRegions,WfmInf-method}
\alias{getSigProbes,WfmInf-method}
\alias{getTwoSided,WfmInf-method}
\alias{getSigGenes,WfmFit,WfmInf-method}
\alias{getNonAnnotatedRegions,WfmFit,WfmInf-method}
\alias{getVarEff,WfmInf-method}
\alias{initialize,WfmInf-method}
\alias{plotWfm,WfmFit,WfmInf-method}
\alias{show,WfmInf-method}

\title{Class "WfmInf"}

\description{class to store outputs from the inference in the wavelet-based transcriptome analysis.}

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
    \item{initialize}{\code{signature(.Object = "WfmInf")}: ... }
    \item{show}{\code{signature(object = "WfmInf")}: ... }
    \item{plotWfm}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
    \item{getSigGenes}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
    \item{getNonAnnotatedRegions}{\code{signature(fit = "WfmFit", inf = "WfmInf")}: ... }
	 }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("WfmInf")
}

\keyword{classes}


