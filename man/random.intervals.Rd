\name{random.intervals}
\alias{random.intervals}
\title{Random intervals}
\usage{
random.intervals(n, M)
}
\arguments{
  \item{n}{the maximal length of an interval}

  \item{M}{a number of intervals}
}
\value{
a matrix with end-points of an interval in each row
}
\description{
Generates endpoints of \code{M} intervals, which are drawn
uniformly with replacements.
}
\examples{
random.intervals(10,100)
}
\seealso{
\code{\link{fixed.intervals}} \code{\link{wbs}}
}

