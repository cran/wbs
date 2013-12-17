\name{wbs}
\alias{wbs}
\alias{wbs.default}
\title{Change-point detection via Wild Binary Segmentation.}
\usage{
wbs(x, ...)

\method{wbs}{default}(x, M = 5000, rand.intervals = TRUE,
  integrated = TRUE, ...)
}
\arguments{
  \item{x}{a vector}

  \item{M}{the number of intervals used in WBS procedure}

  \item{rand.intervals}{a logical variable indicating which
  function is used to generate intervals in the procedure.
  If \code{rand.intervals=TRUE},
  \code{\link{random.intervals}} is applied, otherwise
  \code{\link{fixed.intervals}} is chosen.}

  \item{integrated}{a logical variable indicating the
  version of Wild Binary Segmentation algorithm used. When
  \code{integrated=TRUE}, augmented version of WBS is
  launched, which combines WBS and BS into one.}

  \item{...}{additional arguments that may be passed to wbs
  method.}
}
\value{
an object of class "wbs", which contains the following
fields \item{x}{the vector provided} \item{n}{the length of
\code{x}} \item{M}{the number of intervals used}
\item{rand.intervals}{a logical variable indicating type of
intervals} \item{integrated}{a logical variable indicating
type of WBS procedure} \item{res}{a 5-column matrix with
results, where 's' and 'e' denote start- end points of the
intervals in which change-points candidates 'cpt' have been
found. Column 'CUSUM' contains corresponding value of CUSUM
statistic, 'min.th' is the smallest threshold value for
which change-point candidate would be not added to the set
of estimated change-points.}
}
\description{
Finds change-points for all possible thresholds, using Wild
Binary Segmentation algorithm.
}
\examples{
x <- rnorm(300) + c(rep(1,50),rep(0,250))
w <- wbs(x)
w
plot(w)
cpt <- changepoints(w)
cpt
th <- c(cpt$th,0.7*cpt$th)
cpt <- changepoints(w,th=th)
cpt$cpt.th
}

