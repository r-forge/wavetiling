\name{waveTilingFeatureSet-class}
\Rdversion{1.1}
\docType{class}
\alias{waveTilingFeatureSet-class}
\alias{addPheno,waveTilingFeatureSet-method}
\alias{bgCorrQn,waveTilingFeatureSet-method}
\alias{filterOverlap,waveTilingFeatureSet-method}
\alias{getGroupNames,waveTilingFeatureSet-method}
\alias{getNoGroups,waveTilingFeatureSet-method}
\alias{getReplics,waveTilingFeatureSet-method}
\alias{makeDesign,waveTilingFeatureSet-method}
\alias{selectProbesFromTilingFeatureSet,waveTilingFeatureSet-method}
\alias{wfm.analysis,waveTilingFeatureSet-method}

\title{Class "waveTilingFeatureSet"}

\description{
## %%  ~~ A concise (1-5 lines) description of what the class is. ~~
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("waveTilingFeatureSet", assayData, phenoData, ...)}.
%%  ~~ describe objects here ~~ 
}
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
    \item{addPheno}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{bgCorrQn}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{filterOverlap}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{getGroupNames}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{getNoGroups}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{getReplics}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{makeDesign}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{selectProbesFromTilingFeatureSet}{\code{signature(object = "waveTilingFeatureSet")}: ... }
    \item{wfm.analysis}{\code{signature(object = "waveTilingFeatureSet")}: ... }
	 }
}

\references{
## %%  ~~put references to the literature/web site here~~
}

\author{
## %%  ~~who you are~~
}

\note{
## %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
##
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}

\examples{
showClass("waveTilingFeatureSet")
}

\keyword{classes}
