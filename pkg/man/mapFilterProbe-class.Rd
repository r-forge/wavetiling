\name{mapFilterProbe-class}
\Rdversion{1.1}
\docType{class}
\alias{mapFilterProbe-class}
\alias{getChromosome,mapFilterProbe-method}
\alias{getFilteredIndices,mapFilterProbe-method}
\alias{getPosition,mapFilterProbe-method}
\alias{getStrand,mapFilterProbe-method}
\alias{initialize,mapFilterProbe-method}
\alias{selectProbesFromFilterOverlap,mapFilterProbe-method}
\alias{show,mapFilterProbe-method}

\title{Class "mapFilterProbe"}
\description{
## %%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("mapFilterProbe", filteredIndices, chromosome, position, strand)}.
%%  ~~ describe objects here ~~ 
}
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
    \item{getChromosome}{\code{signature(object = "mapFilterProbe")}: ... }
    \item{getFilteredIndices}{\code{signature(object = "mapFilterProbe")}: ... }
    \item{getPosition}{\code{signature(object = "mapFilterProbe")}: ... }
    \item{getStrand}{\code{signature(object = "mapFilterProbe")}: ... }
    \item{initialize}{\code{signature(.Object = "mapFilterProbe")}: ... }
    \item{selectProbesFromFilterOverlap}{\code{signature(object = "mapFilterProbe")}: ... }
    \item{show}{\code{signature(object = "mapFilterProbe")}: ... }
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
showClass("mapFilterProbe")
}
\keyword{classes}
