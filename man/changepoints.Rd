\name{changepoints}
\alias{changepoints}
\alias{changepoints.sbs}
\alias{changepoints.wbs}
\title{Extract change-points from a sbs or wbs object}
\usage{
changepoints(object, ...)

\method{changepoints}{sbs}(object, th = NULL, th.const = 1.3, Kmax = NULL,
  ...)

\method{changepoints}{wbs}(object, th = NULL, th.const = 1.3, Kmax = 50,
  ssic.penalty = "log", ssic.const = 1.01, ...)
}
\arguments{
  \item{object}{an object of 'wbs' or 'sbs' class}

  \item{th}{a vector of thresholds}

  \item{th.const}{a vector of threshold constants}

  \item{Kmax}{a number of change-points to be detected}

  \item{ssic.penalty}{type of penalty used in sSIC}

  \item{ssic.const}{a constant in penalty in sSIC}

  \item{...}{additional arguments.}
}
\value{
\item{sigma}{Median Absolute Deviation estimate of the
noise level} \item{th}{a vector of thresholds}
\item{no.cpt.th}{the number of change-points detected for
each value of \code{th}} \item{cpt.th}{a list with the
change-points detected for each value of \code{th}}
\item{Kmax}{a maximal number of change-points detected}
\item{ssic.penalty}{type of penalty used in sSIC}
\item{ssic.const}{a constant in penalty in sSIC}
\item{ic.curve}{a list with values of
\code{sSIC},\code{mBIC} and \code{BIC}}
\item{no.cpt.th}{the number of change-points detected for
each information criterion considered} \item{cpt.ic}{a list
with the change-points detected for each information
criterion considered}
}
\description{
The function extracts change-points from \code{object}.
}
\details{
When \code{th=NULL} and \code{Kmax=NULL}, threshold value
is set to \deqn{th=\sigma th.const
\sqrt{2\log(n)}}{th=\sigma*th.const*\sqrt{2\log(n)}}, where
sigma is Median Absolute Deviation estimate of the noise
level.
}

