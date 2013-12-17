#' Print function for 'wbs' objects.
#'
#' @param x an object of class 'wbs'
#' @param ... Ignored
#' @return NULL
#' @method print wbs
#' @S3method print wbs
#' @export
#' @seealso \code{\link{wbs}}

print.wbs <- function(x,...){
	
	cat("Algorithm: ")
	if(x$integrated) cat("Wild Binary Segmentation integreated with standard BS\n")
	else  cat("Wild Binary Segmentation\n")
	
	cat(paste("Number of intervals M =",x$M,"\n"))
	cat("Type of intervals: ")
	if(x$rand.intervals) cat("random\n")
	else  cat("fixed\n")
	cat("Results: \n")
	print(x$res)
}



		
#' Print function for 'sbs' objects.
#'
#' @param x an object of class 'sbs'
#' @param ... Ignored
#' @return NULL
#' @method print sbs
#' @S3method print sbs
#' @export
#' @seealso \code{\link{sbs}}
		
print.sbs <- function(x,...){
	
	cat("Algorithm: standard Binary Segmentation\n")
	cat("Results: \n")
	print(x$res)
}


		