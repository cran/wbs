\docType{package}
\name{wbs-package}
\alias{wbs-package}
\title{Wild Binary Segmentation for multiple change-point detection}
\description{
The package implements Wild Binary Segmentation, a
technique for consistent estimation of the number and
locations of multiple change-points in data. Aditionally, a
fast implementation of standard Binary Segmentation
algorithm is provided.
}
\details{
The main functions of the package are \code{\link{wbs}},
\code{\link{sbs}} and \code{\link{changepoints}}.
}
\examples{
#an example in which standard Binary Segmentation fails to detect change points
x <- rnorm(300)+ c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))

s <- sbs(x)
w <- wbs(x)

cpt.sbs <- changepoints(s)
cpt.sbs

cpt.wbs <- changepoints(w)
cpt.wbs
# in this example, both algorithms work well
x <- rnorm(300) + c(rep(1,50),rep(0,250))

s <- sbs(x)
w <- wbs(x)

cpt.sbs <- changepoints(s)
cpt.sbs

cpt.wbs <- changepoints(w)
cpt.wbs
}
\references{
P. Fryzlewicz (2013), Wild Binary Segmentation for multiple
change-point detection. Under revision.\cr
(\url{http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf})
}

