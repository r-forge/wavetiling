\name{selectProbesFromFilterOverlap}
\alias{selectProbesFromFilterOverlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
select probes (from FilterOverlap)
}

\description{
## %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
selectProbesFromFilterOverlap(object, chromosome, strand=c("forward","reverse"), minPos=min(getPosition(object)), maxPos=max(getPosition(object)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class FilterOverlap
}
  \item{chromosome}{
chromosome
}
  \item{strand}{
strand
}
  \item{minPos}{
minimum position
}
  \item{maxPos}{
maximum Position
}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}

\value{
list
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
\keyword{manip}
