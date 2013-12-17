\name{fixed.intervals}
\alias{fixed.intervals}
\title{Fixed intervals.}
\usage{
fixed.intervals(n, M)
}
\arguments{
  \item{n}{the maximal length of an interval}

  \item{M}{a number of intervals}
}
\value{
a matrix with endpoints of an interval in each row
}
\description{
Generates endpoints of \code{M} intervals.
}
\details{
Function finds the minimal \code{m} such that \code{M} <=
\code{m(m-1)/2}. Then it finds \code{m} approximately
equally-spaced integers and returns all possible intervals
consisting of any two of these points.
}
\examples{
fixed.intervals(10,100)
}
\seealso{
\code{\link{random.intervals}} \code{\link{wbs}}
}

