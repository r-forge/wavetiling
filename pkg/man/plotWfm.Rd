\name{plotWfm}
\alias{plotWfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot Wfm object
}
\description{
plot Wfm fit and detected regions}
\usage{
plotWfm(object, annoFile, minPos, maxPos, trackFeature="exon", overlayFeature=c("gene","transposable_element_gene"), two.strand=TRUE, plotData=TRUE, plotMean=TRUE, tracks=0)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class Wfm
}
  \item{annoFile}{
annotation file
}
  \item{minPos}{
minimum position to plot
}
  \item{maxPos}{
maximum position to plot
}
  \item{trackFeature}{
track feature
}
  \item{overlayFeature}{
overlay feature
}
  \item{two.strand}{
plot one or two strands
}
  \item{plotData}{
plot data or not
}
  \item{plotMean}{
plot mean function or not
}
  \item{tracks}{
number of tracks on plot
}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
nothing returned (check)
}

\references{
## %% ~put references to the literature/web site here ~
}

\author{
Kristof De Beuf <kristof.debeuf@ugent.be>
}

\note{
## %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
## %% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} %check!
