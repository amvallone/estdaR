#' @name theta
#' @rdname theta
#'
#' @title Rank Depcompostion index
#'
#' @description Compute the Regime mobility measure Rey (2004)
#'
#'
#' @param x a \eqn{nxk} matrix with \eqn{l \geq 2} successive columns of a variable are later moments in time
#' @param Regime values corresponding to which regime each observation belongs to
#' @param nsim  number of random spatial permutations to generate for computationally based inference
#'
#' @details For sequence of time periods Theta measures the extent to which rank changes for a variable measured over n locations are in the same direction within mutually exclusive and exhaustive partitions (regimes) of the n locations.
#'
#' @return A data frame
#'
#' @references Rey, S.J. (2004) “Spatial dependence in the evolution of regional income distributions,” in A. Getis, J. Mur and H.Zoeller (eds). Spatial Econometrics and Spatial Statistics. Palgrave, London, pp. 194-213.
#'
#' @examples
#' data(mexico)
#' theta(mexico[,1:7],mexico[,8],999)
#'
#' @export

theta <- function(x,Regime,nsim=NULL){
	rank <- apply(as.matrix(x),2,rank)
	thetas <- rep(0,ncol(x)-1)
	for (j in seq_along(thetas)){
		thetas[j] <- theta.int(rank[,j],rank[,j+1],Regime)
	}
	if(!is.null(nsim)){
		rep <- c()
		for (i in seq_len(nsim)){
			r.thetas <- rep(0,ncol(x)-1)
				for (j in seq_along(r.thetas)){
					r.thetas[j] <- theta.int(rank[,j],rank[,j+1],sample(Regime))
				}
			rep <- rbind(rep,r.thetas)
			}
		larger <- c()
		for (i in seq_along(thetas)){
			larger[i] <- sum(rep[,i]>thetas[i])
		}
		lower <- nsim - larger
		final <- pmin(larger,lower)
	 	p.value<- (final+1)/(nsim+1)
	 	output <- data.frame("theta"=thetas,"p-value"=p.value)
	 	} else {
	 		output <- data.frame("theta"=thetas)
	 	}
	return(output)
}



theta.int <- function(r0,r1,Regime){
	diff <- r1-r0
	class <- unique(Regime)
	aux <- matrix(0,nrow=length(r0),ncol=length(class))
	for (k in seq_along(class)){
		aux[,k] <- Regime == class[k]
	}
	s.col <- colSums(aux*diff)
	out <- sum(abs(s.col))/sum(abs(diff))
	return(out)
}






