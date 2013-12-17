#' @title Change-point detection via standard Binary Segmentation.
#' @description Finds change-points for all possible thresholds, using standard Binary Segmentation algorithm. 
#' @param x a vector
#' @param ... additional arguments that may be passed to sbs method.
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- sbs(x)
#' w
#' plot(w)
#' cpt <- changepoints(w)
#' cpt
#' th <- c(cpt$th,0.7*cpt$th) 
#' cpt <- changepoints(w,th=th)
#' cpt
#' @rdname sbs
#' @export
#' @return an object of class "sbs", which contains the following fields
#' \item{x}{the vector provided}
#' \item{n}{the length of \code{x}}
#' \item{res}{a 5-column matrix with results, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'cpt' have been found.
#' Column 'CUSUM' contains corresponding value of CUSUM statistic, 'min.th' is the smallest 
#' threshold value for which change-point candidate would be not added to the set of estimated 
#' change-points.}

sbs <- function(x, ...)  UseMethod("sbs")

#' @method sbs default
#' @S3method sbs default
#' @rdname sbs

sbs.default <- function(x, ...){
	results <- list()	
	results$x <- as.numeric(x)
	results$n <- length(results$x)

	
	if(results$n <2) stop("x should contain at least two elements")
	if(NA%in%results$x) stop("x vector cannot contain NA's")
	if(var(x)==0) stop("x is a constant vector, change-point detection is not needed")
	tmp <- .C("bs_rec_wrapper", as.double(results$x), as.integer(results$n), as.double(rep(0,5*(results$n-1))))[[3]]
	results$res <- matrix(tmp,ncol=5)		
	
	colnames(results$res) <- c("s","e","cpt","CUSUM","min.th")
	
	class(results) <- "sbs"
	return(results)
	
}
