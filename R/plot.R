#' @title Plot function for 'sbs' object.
#' @description Plots data with fitted piecewiese constant function, equal to the mean
#' between specified change-points.
#' @details When \code{cpt} is omitted, the function automatically finds change-points 
#' using \code{changepoints} function.
#' @method plot sbs
#' @S3method plot sbs
#' @export 
#' @param x an object of class 'sbs', returned by \code{\link{sbs}}
#' @param cpt a vector with localisations of change-points
#' @param ... other parameters which may be passed to \code{plot} and \code{changepoints}
#' @seealso \code{\link{sbs}}  \code{\link{changepoints}}

plot.sbs <- function(x,cpt,...){
	ts.plot(x$x,ylab="x",...)
	
	if(missing(cpt)){
		w.cpt <- changepoints(x,...)
		print
		means <- mean.from.cpt(x$x,w.cpt$cpt.th[[1]])
	}else{
		means <- mean.from.cpt(x$x,cpt)
	}
	
	lines(x=means,type="l",col="red")
	title("Fitted piecewise constant function")
	
}

#' @title Plot function for 'wbs' object.
#' @description Plots data with fitted piecewiese constant function, equal to the mean
#' between specified change-points.
#' @details When \code{cpt} is omitted, the function automatically finds change-points 
#' using \code{changepoints} function.
#' @method plot wbs
#' @S3method plot wbs
#' @export 
#' @param x an object of class 'wbs', returned by \code{\link{wbs}}
#' @param cpt a vector with localisations of change-points
#' @param ... other parameters which may be passed to \code{plot} and \code{changepoints}
#' @seealso \code{\link{wbs}}  \code{\link{changepoints}}


plot.wbs <- function(x,cpt,...){
	
	if(missing(cpt)) plot.sbs(x,...)
	else plot.sbs(x,cpt,...)
	
}

