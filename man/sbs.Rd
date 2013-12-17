\name{sbs}
\alias{sbs}
\alias{sbs.default}
\title{Change-point detection via standard Binary Segmentation.}
\usage{
sbs(x, ...)

\method{sbs}{default}(x, ...)
}
\arguments{
  \item{x}{a vector}

  \item{...}{additional arguments that may be passed to sbs
  method.}
}
\value{
an object of class "sbs", which contains the following
fields \item{x}{the vector provided} \item{n}{the length of
\code{x}} \item{res}{a 5-column matrix with results, where
's' and 'e' denote start- end points of the intervals in
which change-points candidates 'cpt' have been found.
Column 'CUSUM' contains corresponding value of CUSUM
statistic, 'min.th' is the smallest threshold value for
which change-point candidate would be not added to the set
of estimated change-points.}
}
\description{
Finds change-points for all possible thresholds, using
standard Binary Segmentation algorithm.
}
\examples{
x <- rnorm(300) + c(rep(1,50),rep(0,250))
w <- sbs(x)
w
plot(w)
cpt <- changepoints(w)
cpt
th <- c(cpt$th,0.7*cpt$th)
cpt <- changepoints(w,th=th)
cpt
}

