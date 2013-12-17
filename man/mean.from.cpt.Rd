\name{mean.from.cpt}
\alias{mean.from.cpt}
\title{Evaluates the mean between cpt.}
\usage{
\method{mean}{from.cpt}(x, cpt = NULL, ...)
}
\arguments{
  \item{x}{a vector}

  \item{cpt}{a vector localisations of change-points}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
returns a vector of the same lenght as \code{x}, piecewise
constant and equal to the mean between change-points given
in \code{cpt}
}
\description{
Finds the mean between given cpt.
}
\examples{
x <- rnorm(100)+c(rep(-1,50),rep(1,50))
cpt <- 50
mean.from.cpt(x,cpt)
w <- wbs(x)
cpt <- changepoints(w)
mean.from.cpt(x,cpt=cpt$cpt.ic$sbic)
}

