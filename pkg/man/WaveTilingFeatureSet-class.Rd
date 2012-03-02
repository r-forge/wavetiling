\name{WaveTilingFeatureSet-class}
\Rdversion{1.1}
\docType{class}
\alias{WaveTilingFeatureSet-class}
\alias{initialize,WaveTilingFeatureSet-method}
\alias{show,WaveTilingFeatureSet-method}
%\alias{pmStrand,WaveTilingFeatureSet-method}
\alias{addPheno,WaveTilingFeatureSet-method}
\alias{bgCorrQn,WaveTilingFeatureSet-method}
\alias{filterOverlap,WaveTilingFeatureSet-method}
\alias{getGroupNames,WaveTilingFeatureSet-method}
\alias{getNoGroups,WaveTilingFeatureSet-method}
\alias{getReplics,WaveTilingFeatureSet-method}
\alias{selectProbesFromTilingFeatureSet,WaveTilingFeatureSet-method}
\alias{wfm.fit,WaveTilingFeatureSet-method}

\title{Class "WaveTilingFeatureSet"}

\description{Class to store expression and phenotypic data from a tiling array experiment, used as input for the wavelet-based transcriptome analysis.}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("WaveTilingFeatureSet")}.}

\section{Slots}{
  \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"} ~~ }
    \item{\code{intensityFile}:}{Object of class \code{"character"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{TilingFeatureSet}"}, directly.
Class \code{"\linkS4class{FeatureSet}"}, by class "TilingFeatureSet", distance 2.
Class \code{"\linkS4class{NChannelSet}"}, by class "TilingFeatureSet", distance 3.
Class \code{"\linkS4class{eSet}"}, by class "TilingFeatureSet", distance 4.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "TilingFeatureSet", distance 5.
Class \code{"\linkS4class{Versioned}"}, by class "TilingFeatureSet", distance 6.
}

\section{Methods}{
  \describe{
    \item{addPheno}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{bgCorrQn}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{filterOverlap}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{getGroupNames}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{getNoGroups}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{getReplics}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{selectProbesFromTilingFeatureSet}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
    \item{wfm.fit}{\code{signature(object = "WaveTilingFeatureSet")}: ... }
	 }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("WaveTilingFeatureSet")
}

\keyword{classes}
