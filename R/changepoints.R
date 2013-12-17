#' @title Extract change-points from a sbs or wbs object
#' @description The function extracts change-points from \code{object}.
#' @details When \code{th=NULL} and \code{Kmax=NULL},  threshold value is set to
#' \deqn{th=\sigma th.const \sqrt{2\log(n)}}{th=\sigma*th.const*\sqrt{2\log(n)}},
#' where sigma is Median Absolute Deviation estimate of the noise level.
#' @param object an object of 'wbs' or 'sbs' class
#' @param th a vector of thresholds
#' @param th.const a vector of threshold constants
#' @param Kmax a number of change-points to be detected 
#' @param ssic.penalty type of penalty used in sSIC
#' @param ssic.const a constant in penalty in sSIC 
#' @param ... additional arguments.
#' @return  
#' \item{sigma}{Median Absolute Deviation estimate of the noise level}
#' \item{th}{a vector of thresholds}
#' \item{no.cpt.th}{the number of change-points detected for each value of \code{th}}
#' \item{cpt.th}{a list with the change-points detected for each value of \code{th}}
#' \item{Kmax}{a maximal number of change-points detected}
#' \item{ssic.penalty}{type of penalty used in sSIC}
#' \item{ssic.const}{a constant in penalty in sSIC}
#' \item{ic.curve}{a list with values of \code{sSIC},\code{mBIC} and \code{BIC}}
#' \item{no.cpt.th}{the number of change-points detected for each information criterion considered}
#' \item{cpt.ic}{a list with the change-points detected for each information criterion considered}
#' @rdname changepoints
#' @export

changepoints  <- function(object,...) UseMethod("changepoints")


#' @method changepoints sbs
#' @S3method changepoints sbs
#' @rdname changepoints

changepoints.sbs <- function(object,th=NULL,th.const=1.3,Kmax=NULL,...){
	
	th <- as.numeric(th)
	th.const <- as.numeric(th.const)
	Kmax <- as.integer(Kmax)
	
	w.cpt <- list()
	
	w.cpt$sigma <- mad(diff(object$x)/sqrt(2))
	
	if(length(th)) w.cpt$th <- as.numeric(th)
	else if(length(Kmax)){
		if(Kmax[1] < 1) stop("Kmax should be at least 1")
		Kmax <- as.integer(min(nrow(object$res),Kmax))
		if(Kmax < nrow(object$res)){
			tmp <- sort(object$res[,5],decreasing=TRUE)
			w.cpt$th <- tmp[Kmax]-min(diff(sort(unique(tmp),decreasing=FALSE)))/2
		}else w.cpt$th <- 0;
		
	}else{
		w.cpt$th <-  th.const * w.cpt$sigma *sqrt(2*log(object$n))
	}
	

	
	th.l <- length(w.cpt$th)
	
	if(th.l) w.cpt$th <- sort(w.cpt$th)
	else stop("th vector should contain at least one value")
	
	if(NA%in%w.cpt$th) stop("th cannot contain NA's")
	if(NA%in%w.cpt$Kmax) stop("Kmax cannot be NA")
	
	if(length(Kmax)>1){
		Kmax <- Kmax[1]
		warning("Kmax could contain only one number")
	}
	

	
	res.tmp <- matrix(object$res[order(abs(object$res[,4]),decreasing=TRUE),c(3,5)],ncol=2)
	
	w.cpt$no.cpt.th <- as.integer(rep(0,th.l))
	
	for(i in 1:th.l)
		if(nrow(res.tmp)){
			
			res.tmp <- matrix(res.tmp[res.tmp[,2]>w.cpt$th[i],],ncol=2)
			if(nrow(res.tmp)){
				w.cpt$cpt.th[[i]] <- res.tmp[,1]
				w.cpt$no.cpt.th[i] <- length(w.cpt$cpt.th[[i]])
			}else w.cpt$cpt.th[[i]] <- NA
		}else w.cpt$cpt.th[[i]] <- NA
		
		
	w.cpt$Kmax <- as.integer(max(w.cpt$no.cpt.th))	
	
	
		
	return(w.cpt)
	
}

#' @method changepoints wbs
#' @S3method changepoints wbs
#' @rdname changepoints

changepoints.wbs <- function(object,th=NULL,th.const=1.3, Kmax=50, ssic.penalty="log", ssic.const=1.01,...){
	

	w.cpt <- changepoints.sbs(object,th,th.const=th.const,Kmax=NULL)
	w.cpt$ssic.penalty <- as.character(ssic.penalty)
	w.cpt$ssic.const <- as.numeric(ssic.const)

	
	
	#ic part 
	if(Kmax==0 || object$n<=3) {
		w.cpt$cpt.ic$ssic <- w.cpt$cpt.ic$bic <- 	w.cpt$cpt.ic$mbic <- NA
		w.cpt$no.cpt.ic[c("ssic","mbic","bic")] <- as.integer(c(0,0,0))
		w.cpt$ic.curve <- list()
		w.cpt$ic.curve$bic<-w.cpt$ic.curve$ssic<-w.cpt$ic.curve$mbic <- NA
	}else{
		
		w.cpt$Kmax <- as.integer(Kmax)
		cpt.cand <- changepoints.sbs(object,th,th.const=1.3,Kmax=w.cpt$Kmax)$cpt.th[[1]]		
		if(NA%in%cpt.cand) stop("no change-poinst found, choose larger Kmax")
		
		if(length(cpt.cand)>min(Kmax,object$n-2)) cpt.cand <- cpt.cand[1:min(Kmax,object$n-2)]
		len.cpt<- length(cpt.cand)
		

		
		w.cpt$ic.curve$ssic <- w.cpt$ic.curve$mbic <- w.cpt$ic.curve$mbic <- rep(0,len.cpt+1) 
		
		if (w.cpt$ssic.penalty == "log") pen <- log(object$n)^w.cpt$ssic.const
		else if (w.cpt$ssic.penalty == "power") pen <- object$n^w.cpt$ssic.const
		else{
			warning("unknown ssic penalty specified, setting \"log\"")
			w.cpt$ssic.penalty = "log"
		} 

		
		if(len.cpt) for(i in len.cpt:1){
			means <- mean.from.cpt(object$x,cpt.cand[1:i])
			min.log.lik <- object$n/2 * log(sum((object$x - means)^2)/object$n)
			w.cpt$ic.curve$ssic[i+1] <- min.log.lik + i * pen
			w.cpt$ic.curve$bic[i+1] <- min.log.lik + i * log(object$n)
			w.cpt$ic.curve$mbic[i+1] <-min.log.lik + 3/2 * i * log(object$n) + 1/2 * sum(log(diff(c(0,sort(cpt.cand[1:i]),object$n))/object$n))
		}
		
		w.cpt$ic.curve$ssic[1]<- w.cpt$ic.curve$mbic[1] <- w.cpt$ic.curve$bic[1] <- min.log.lik<- object$n/2 * log(var(object$x))
		
		
		criteria <- c("bic","ssic","mbic")
		w.cpt$cpt.ic <- list()
		
		for(i in 1:3){
			tmp <- quantile(which.min(w.cpt$ic.curve[[criteria[i]]]), .5, type=3)
			
			if(tmp ==1){
				w.cpt$cpt.ic[[criteria[i]]] <- NA
				w.cpt$no.cpt.ic[criteria[i]] <- as.integer(0)
			}else{
				w.cpt$cpt.ic[[criteria[i]]] <- cpt.cand[1:(tmp-1)]
				w.cpt$no.cpt.ic[criteria[i]] <- as.integer(tmp-1);
			}
		}
		
		
		
	}			
	
	return(w.cpt)
	
}