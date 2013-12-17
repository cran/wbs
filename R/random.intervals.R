#' Random intervals
#'
#' Generates endpoints of \code{M} intervals, which are drawn uniformly with replacements.
#' 
#' @param n the maximal length of an interval
#' @param M a number of intervals 
#' @return a matrix with end-points of an interval in each row
#' @examples
#' random.intervals(10,100)
#' @export random.intervals
#' @seealso \code{\link{fixed.intervals}} \code{\link{wbs}}

random.intervals <-
		function(n,M) {
	
	n <- as.integer(n)
	M <- as.integer(M)
	intervals <- matrix(0,nrow=M,ncol=2)
	for (i in 1:M) {
		ind <- floor(runif(2) * (n-1)) + c(1,2)
		s <- min(ind)
		e <- max(ind)
		if (e == s) e <- e+1
		intervals[i,1:2] <- c(s,e)

	}
	intervals
}


#random.intervals <-
#		function(n,M) {
#	
#	n <- as.integer(n)
#	M <- as.integer(M)
#	intervals <- matrix(0,nrow=M,ncol=2)
#	intervals[,1] <- ceiling(runif(M)*(n-1))
#	intervals[,2] <- intervals[,1]+ ceiling(runif(M)*(n-intervals[,1]))
#	
#}



