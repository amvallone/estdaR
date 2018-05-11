#' @name sp.mkv
#' @rdname sp.mkv
#'
#' @title Spatial Markov transition probability matrix
#' @description Compute the Spatial Markov transition probability matrix (Rey,2001)
#'
#' @param m numerical matrix of n spatial unit ans t time periods
#' @param W an objet of listw class.
#' @param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
#' @param fixed logical, if it is TRUE the data are pooled over space and time and the quintiles calculated for the pooled data
#'
#' @details for later..
#' @return a list contaning the markov's trantiiona matrix and the markov's transition probability matrix
#'
#' @references Rey, S.J. (2001) “Spatial empirics for economic growth and convergence”, 34 Geographical Analysis, 33, 195-214.
#'
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' pci <- data[,10:90]
#' rpci <- pci/matrix(1,dim(pci))%*%colMeans(pci)
#' w1queen <- nb2listw(poly2nb(us48))
#' sp.mkv(rpci,w1queen)
#'
#' @export


sp.mkv <- function(m,W,classes=5,fixed=TRUE){
	lag<-matrix(0,nrow(m),ncol(m))
	for(j in 1:ncol(m)){
		lag[,j]<-lag.listw(W,m[,j])
	}
	if(fixed==TRUE){
		auxl <- unlist(c(lag))
		auxx <- unlist(c(m))
		auxl <- discret(auxl,classes=classes)
		auxx <- discret(auxx,classes=classes)
		lx <- matrix(auxl,dim(lag))
		x <- matrix(auxx,dim(m))
	} else {
		lx <- apply(lag,2,discret,classes=classes)
		x <- apply(m,2,discret,classes=classes)
	}
	class<-as.numeric(unique(as.factor(lx)))
    class<-sort(class)
    mm <- array(0,dim = c(length(class),length(class),length(class)),dimnames=list(class,class,paste("Lag",class)))
	for (j in 1:ncol(x)-1){
		for(i in 1:nrow(x)){
        		mm[x[i,j], x[i,j+1],lx[i,j]] <- mm[x[i,j], x[i,j+1],lx[i,j]] + 1
      		}
    }
    mm.p <- array(0,dim = c(length(class),length(class),length(class)),dimnames=list(class,class,paste("Lag",class)))
    for(k in seq_along(class)){
    		mm.p[,,k] <- mm[,,k]/rowSums(mm[,,k])
    }
  output <- list("Probabilities" = mm.p, "Transitions" = mm)
  return(output)
 }
