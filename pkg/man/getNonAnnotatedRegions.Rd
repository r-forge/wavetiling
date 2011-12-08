\name{getNonAnnotatedRegions}
\alias{getNonAnnotatedRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get non-annotated regions
}

\description{
## %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
getNonAnnotatedRegions(object,annoFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class Wfm
}
  \item{annoFile}{
annotation file
}
}

\details{
## %%  ~~ If necessary, more details than the description above ~~
}
\value{
list of GRanges object with the non-annotated regions, both for no annotation overlap on the strand under consideration, and for no annotation overlap on both strands
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
