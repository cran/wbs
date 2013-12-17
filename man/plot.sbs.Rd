\name{plot.sbs}
\alias{plot.sbs}
\title{Plot function for 'sbs' object.}
\usage{
\method{plot}{sbs}(x, cpt, ...)
}
\arguments{
  \item{x}{an object of class 'sbs', returned by
  \code{\link{sbs}}}

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
\code{\link{sbs}} \code{\link{changepoints}}
}

