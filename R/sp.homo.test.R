#' @name sp.homo.test
#' @rdname sp.homo.test
#'
#' @title Test for homogeneity of Markov transition probabilities across regimes.
#' @description Performs the homogenity across space test for spatial markov trasntion matrix basis on Rey et al, (2016)
#'
#' @param x numerical matrix of n spatial unit ans t time periods
#' @param W an objet of listw class.
#' @param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
#' @param fixed logical, if it is TRUE the data are pooled over space and time and the quintiles calculated for the pooled data
#'
#' @details For later...
#'
#' @return A list coantaning the Q statistic, the LR statistic and matrix use as null hypotesis in the test.
#'
#' @references S. J. Rey, W. Kang, and L. Wolf (2016) “The properties of tests for spatial effects in discrete Markov chain models of regional income distribution dynamics,” Journal of  Geographical Systems, vol. 18, no. 4, pp. 377–398.
#'
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' pci <- data[,10:90]
#' rpci <- pci/matrix(1,dim(pci))%*%colMeans(pci)
#' w1queen <- nb2listw(poly2nb(us48))
#' sp.homo.test(rpci,w1queen)
#'
#' @export

sp.homo.test <- function(x,W,classes=5,fixed=TRUE){
	n <- classes
	sm <- sp.mkv(x,W,classes=classes,fixed=fixed)
	M <- sm[[2]]
	B <- matrix(0,n,n)
	T.M <- apply(M,c(1,2),sum)
	tot <- sum(T.M)
	n_i <- rowSums(T.M)
	A_i <- rowSums(T.M>0)
	A_im <- matrix(0,n,n)
	p_ij <- diag(1/(n_i+sum(n_i==0)))%*%T.M
	den <- p_ij+1*(p_ij==0)
	b_i <- A_i*0
	p_ijm <- M*0
	Q <- 0
	L.R <- 0
	q.table <- M*0
	LR <- M*0
	k<-1
	for (j in 1:n){
		m <- M[,,j]
		nim<-rowSums(m)
		B[,k] <- 1*(nim>0)
		b_i <- b_i + 1 * (nim > 0)
		p_ijm[,,k] <- diag(1/(nim+sum(nim==0))) %*% m
		num <-  (p_ijm[,,k]-p_ij)^2
		ratio <-  num / den
		qijm <- diag(nim) %*% ratio
		q.table[,,k] <- qijm
		Q <- Q + sum(qijm)
		mask = (m > 0) * (p_ij > 0)
		A_im[, k] = rowSums(m > 0)
		unmask <- 1 * (mask == 0)
		lr.ratio <- (mask * p_ijm[,,k] + unmask) / (mask * p_ij + unmask)
		lr = m * log(lr.ratio)
		L.R <- L.R +sum(lr)
		LR[,,k]<-2 * lr
		k+1
	}
	dof <- as.integer(sum((b_i - 1) * (A_i - 1)))
	L.R <- L.R * 2
	pQ <- 1-pchisq(Q,dof)
	pL.R <- 1-pchisq(L.R,dof)
	testQ <- data.frame("Q"=round(Q,2),"p-value"=round(pQ,4),"d.o.f" = round(dof,0))
	testLR <- data.frame("LR"=round(L.R,2),"p-value"=round(pL.R,4),"d.o.f"= round(dof,0))
	out <- list("Q"=testQ, "LR"=testLR, "NULL"=round(p_ij,4))
	return(out)
}

