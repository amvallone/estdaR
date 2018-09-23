#' @name sp.tau
#' @rdname sp.tau
#'
#' @title Global Indicators of Mobility Association (GIMA)
#' @description Compute the Global Indicators of Mobility Association index Rey (2016)
#'
#' @param x rank variable
#' @param y rank variable.
#' @param W an object of class \code{listw}
#' @param perm  number of random spatial permutations for calculation of pseudo p-values, the default value is NULL.
#'
#' @details The Global Indicators of Mobility Association (GIMA) is based on the Spatial Tau indicator (Rey,2004). The implementation is a two step algorith based on the Rey's implementation (Rey,2014)
#'
#' @return a vector coantaning \itemize{
#'    \item The GIMA statistic value
#'    \item The GIMA's pseudo p value, only available for perm >0
#'	  \item The number of spatial Concordant pairs
#'	  \item The number of spatial Discordant pairs
#'    \item The Tau statistic value
#'    \item The tau's p value
#'	  \item The number of Concordant pairs
#' 	  \item The number of Discordant pairs
#'	  \item The number of extra x pairs. An extra x pair is a pair which \eqn{sgn(x_i - x_j)=0}
#'	  \item The number of extra y pairs. An extra y pair is a pair which \eqn{sgn(y_i - y_j)=0}
#'    }
#' @examples
#' data(mexico)
#' n <- nrow(mexico)
#' w <- matrix(0,n,n)
#' for(i in 1:n){
#'		w[,i] <- 1*(mexico$Regime[i]==mexico$Regime)
#' }
#' diag(w) <- 0
#' W <- mat2listw(w,style="B")
#' a <- lapply(1:(ncol(mexico)-2),function(i) sp.tau(mexico[,i],mexico[,(i+1)],W,999))
#' b<-do.call(rbind.data.frame,a)
#' colnames(b) <- names(a[[1]])
#' b <- round(b,3)
#'
#'
#' @export




sp.tau <- function(x,y,W,perm=NULL){
	f.step <- tau(x,y)
	Tau <- f.step["Tau"]
	Tau_p <- f.step["pval"]
	concordant <- f.step["Concordant"]
	discordant <- f.step["Discordant"]
	extraX <- f.step["ExtraX"]
	extraY <- f.step["ExtraY"]
	res <- int.tau(x,y,W)
	if(!is.null(perm)){
		taus <- rep(0,perm)
		ids <- seq_along(x)
		for (r in 1:perm){
			rids <- sample(ids)
			taus[r] <- int.tau(x[rids],y[rids],W)[1]
		}
		above <- taus >= res[1]
		larger <- sum(above)
		psim <- (larger + 1)/ (perm + 1)
		if (psim > 0.5) psim <- (perm - larger +1) / (perm + 1)
		out <- c("sp.tau"=res["tau_g"],"pval"=psim,"sp.Concordant"=res["gc"],"sp.Discordant"=res["gd"],f.step)
	} else {
		c("sp.tau"=res["tau_g"],"sp.Concordant"=res["gc"],"sp.Discordant"=res["gd"],f.step)
	}
	out
}

int.tau<-function(x,y,W){
	n1 <- n2 <- iS <- gc <- 0
	for(i in seq_along(W$neighbours)){
		neg <- unlist(W$neighbours[i])
		xi <- x[i]
		yi <- y[i]
		for(k in seq_along(neg)){
			j <- neg[k]			
			if(i < j){
				xj <- x[j]
				yj <- y[j]
				dx <- xi - xj
				dy <- yi - yj
				dxdy <- dx * dy
				if(dxdy!=0){
					n1 <- n1 + 1
					n2 <- n2 + 1
					if(dxdy > 0){
						gc <- gc + 1
						iS <- iS + 1
					} else {
						iS <- iS - 1
					} 
				} else {
					if (dx !=0) n1 <- n1 +1
					if (dy !=0) n2 <- n2 +1
				}
			}
		}
	}
	tau_g <- iS /(sqrt(n1) * sqrt(n2))
	gd <- gc - iS
	out <- c("tau_g"=tau_g,"gc"=gc,"gd"=gd)
	return(out)
}


