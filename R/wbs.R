#' @title Wild Binary Segmentation for multiple change-point detection
#' @description The package implements Wild Binary Segmentation, a technique for consistent estimation of
#' the number and locations of multiple change-points in data. Aditionally, a fast implementation of standard Binary Segmentation
#' algorithm is provided.
#' @details The main functions of the package are \code{\link{wbs}}, \code{\link{sbs}} and \code{\link{changepoints}}.
#' @references P. Fryzlewicz (2013), Wild Binary Segmentation for multiple change-point detection. Under revision.\cr 
#'  (\url{http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf})
#' @docType package
#' @useDynLib wbs
#' @name wbs-package
#' @examples
#' #an example in which standard Binary Segmentation fails to detect change points
#' x <- rnorm(300)+ c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' 
#' s <- sbs(x)
#' w <- wbs(x)
#' 
#' cpt.sbs <- changepoints(s)
#' cpt.sbs
#' 
#' cpt.wbs <- changepoints(w)
#' cpt.wbs
#' # in this example, both algorithms work well
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' 
#' s <- sbs(x)
#' w <- wbs(x)
#' 
#' cpt.sbs <- changepoints(s)
#' cpt.sbs
#' 
#' cpt.wbs <- changepoints(w)
#' cpt.wbs


NULL

#' @title Change-point detection via Wild Binary Segmentation.
#' @description Finds change-points for all possible thresholds, using Wild Binary Segmentation algorithm. 
#' @param x a vector
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' w
#' plot(w)
#' cpt <- changepoints(w)
#' cpt
#' th <- c(cpt$th,0.7*cpt$th) 
#' cpt <- changepoints(w,th=th)
#' cpt$cpt.th
#' @rdname wbs
#' @export

wbs <- function(x, ...)  UseMethod("wbs")

#' @method wbs default
#' @S3method wbs default
#' @rdname wbs
#' @param M the number of intervals used in WBS procedure
#' @param rand.intervals a logical variable indicating which function is used to generate intervals in the procedure.
#' If \code{rand.intervals=TRUE}, \code{\link{random.intervals}} is applied, otherwise \code{\link{fixed.intervals}} is chosen.      
#' @param integrated a logical variable indicating the version of Wild Binary Segmentation algorithm used. When \code{integrated=TRUE}, 
#' augmented version of WBS is launched, which combines WBS and BS into one. 
#' @param ... additional arguments that may be passed to wbs method.
#' @return an object of class "wbs", which contains the following fields
#' \item{x}{the vector provided}
#' \item{n}{the length of \code{x}}
#' \item{M}{the number of intervals used}
#' \item{rand.intervals}{a logical variable indicating type of intervals}
#' \item{integrated}{a logical variable indicating type of WBS procedure}
#' \item{res}{a 5-column matrix with results, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'cpt' have been found.
#' Column 'CUSUM' contains corresponding value of CUSUM statistic, 'min.th' is the smallest 
#' threshold value for which change-point candidate would be not added to the set of estimated 
#' change-points.}

 

wbs.default <- function(x, M=5000,  rand.intervals = TRUE,integrated=TRUE,...){
	results <- list()	
	results$x <- as.numeric(x)
	results$n <- length(results$x)
	results$M <- as.integer(M)
	results$integrated <- as.logical(integrated)
	results$rand.intervals <- as.logical(rand.intervals)
	results$res <- matrix(nrow=0,ncol=5)
	
	if(results$n <2) stop("x should contain at least two elements")
	if(NA%in%results$x) stop("x vector cannot contain NA's")
	if(var(x)==0) stop("x is a constant vector, change-point detection is not needed")
	
	if(is.na(results$M)) stop("M cannot be NA")
	if(length(results$M)> 1)  stop("M should be a single integer")
	if(results$M<0)  stop("M should be an integer > 0")
	
	if(results$rand.intervals) intervals <-  matrix(random.intervals(results$n,results$M),ncol=2)
	else intervals <- matrix(fixed.intervals(results$n,results$M),ncol=2)
	
	if(results$integrated){
		tmp <- .C("wbs_int_rec_wrapper", as.double(results$x), as.integer(results$n), as.double(rep(0,5*(results$n-1))), as.integer(intervals), as.integer(results$M))[[3]]
		results$res <- matrix(tmp,ncol=5)	
	}else{
		tmp <- .C("wbs", as.double(results$x), as.integer(results$n), as.double(rep(0,5*results$M)), as.integer(intervals), as.integer(results$M))[[3]]
		
		sorted <- matrix(tmp,ncol=5)
		sorted <- matrix(sorted[order(sorted[,5],decreasing=TRUE),],ncol=5)
		
		results$res <- matrix(0, nrow=results$n-1,ncol=5)
		
		j <- 1
		
		while (nrow(sorted)) {
			results$res[j,] <- sorted[1,]
			cpt.cand <- sorted[1,3]
			sorted <- matrix(sorted[!(sorted[,1] <= cpt.cand & sorted[,2] > cpt.cand),], ncol=5)
			j <- j + 1
		}	
		
		results$res <- matrix(results$res[1:(j-1),],ncol=5)
		
	}
	
	
	colnames(results$res) <- c("s","e","cpt","CUSUM","min.th")
	
	class(results) <- "wbs"
	return(results)
	
}


