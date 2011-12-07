\name{genomeInfo-class}
\Rdversion{1.1}
\docType{class}
\alias{genomeInfo-class}
\alias{initialize,genomeInfo-method}

\title{Class "genomeInfo"}
\description{
## %%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("genomeInfo", chromosome, strand, minPos, maxPos)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{chromosome}:}{Object of class \code{"vector"} ~~ }
    \item{\code{strand}:}{Object of class \code{"character"} ~~ }
    \item{\code{minPos}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{maxPos}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "genomeInfo")}: ... }
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
##%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("genomeInfo")
}
\keyword{classes}
