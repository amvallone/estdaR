#' @name join.d
#' @rdname join.d
#'
#' @title Test the independence in the dynamics of a variable and its neighbors
#' @description Performe a Chi square test for test that dynamics of a variable is independent of dynamics of its neighbors.
#'
#' @param x numerical matrix of n spatial unit ans t time periods
#' @param W an objet of listw class.
#'
#' @details The test decompose the LISA Markov chain in a pair of chains, one for the city and other for the neighbors, each chain has two states H and L and under the null of independence test the co-movement of the chains.
#'
#' @return a list coantaning \itemize{
#'    \item the Chi square statistcas and its p value
#'    \item The LISA markov transition matrix ander the null
#'    }
#'
#' @examples
#' data(us48)
#' join.d(data[,9:48],w1rook)
#'
#' @export

join.d <- function(x,W){
	t<-dim(x)[2L]
  n<-dim(x)[1L]
	x.mean <- colMeans(x)
	x.bar <- x/x.mean
	x.L<-matrix(0,nrow=n,ncol=t)
		for (i in 1:t){
			x.L[,i]<-lag.listw(W,x[,i])
		}
	x.L.mean <- colMeans(x.L)
	x.L.bar <- x.L/(x.L.mean)
	h.x <- apply(x.bar > 1,2,as.numeric) + 1 # 2 are a h and 1 a low
	h.x.l <- apply(x.L.bar > 1,2,as.numeric) + 1
	#################################################
	m.x <- mkv.int(h.x)
	m.x <- m.x/rowSums(m.x)
	m.x.l <- mkv.int(h.x.l)
	m.x.l <- m.x.l/rowSums(m.x.l)
	A <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0),4,4,byrow=TRUE)
	kp <- m.x %x% m.x.l
	aux <- A %*% kp %*% t(A)
	trans <- lisamkv(x,W)$lisamatrix
	t.trans <- rowSums(trans)
	t.h0 <- diag(t.trans) %*% aux
	chi <- m.chi(trans,t.h0)
	output <- list("Chi2"=chi,"Expected"=t.h0)
	return(output)
}

