% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penalties.R
\name{bic.penalty}
\alias{bic.penalty}
\title{Bayesian Information Criterion penalty term}
\usage{
bic.penalty(n, cpt)
}
\arguments{
\item{n}{the number of observations}

\item{cpt}{a vector with localisations of change-points}
}
\value{
the penalty term \eqn{k\log(n)}{k * log(n)} where  \eqn{k}{k} denotes the number of elements in \code{cpt}
}
\description{
The function evaluates the penalty term for the standard Bayesian Information Criterion applied to the change-point detection problem. This routine is typically not called directly by the user; its name can be passed as an argument to \code{\link{changepoints}}.
}
\examples{
x <- rnorm(300) + c(rep(1,50),rep(0,250))
w <- wbs(x)
w.cpt <- changepoints(w,penalty="bic.penalty")
w.cpt$cpt.ic
x <- rnorm(300) + c(rep(1,50),rep(0,250))
w <- wbs(x)
w.cpt <- changepoints(w,penalty="bic.penalty")
w.cpt$cpt.ic
}
