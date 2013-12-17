\name{plot.wbs}
\alias{plot.wbs}
\title{Plot function for 'wbs' object.}
\usage{
\method{plot}{wbs}(x, cpt, ...)
}
\arguments{
  \item{x}{an object of class 'wbs', returned by
  \code{\link{wbs}}}

  \item{cpt}{a vector with localisations of change-points}

  \item{...}{other parameters which may be passed to
  \code{plot} and \code{changepoints}}
}
\description{
Plots data with fitted piecewiese constant function, equal
to the mean between specified change-points.
}
\details{
When \code{cpt} is omitted, the function automatically
finds change-points using \code{changepoints} function.
}
\seealso{
\code{\link{wbs}} \code{\link{changepoints}}
}

