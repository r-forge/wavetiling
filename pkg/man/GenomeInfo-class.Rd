\name{GenomeInfo-class}
\Rdversion{1.1}
\docType{class}
\alias{GenomeInfo-class}
\alias{initialize,GenomeInfo-method}
\alias{show,GenomeInfo-method}

\title{Class "genomeInfo"}
\description{class to store the genomic info from a \code{WfmFit}-class object}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("GenomeInfo", chromosome, strand, minPos, maxPos)}.}

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
    \item{initialize}{\code{signature(.Object = "GenomeInfo")}: ... }
    \item{show}{\code{signature(.Object = "GenomeInfo")}: ... }
	 }
}

\author{Kristof De Beuf <kristof.debeuf@ugent.be>}

\examples{
showClass("GenomeInfo")
}

\keyword{classes}
