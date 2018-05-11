#' @name shorrock
#' @rdname shorrock
#'
#' @title Shorrock's mobility measure
#' @description Compute the Shorrock's mobility measure
#'
#' @param x Markov probability transition matrix.
#'
#' @return a vector
#'
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' pci <- data[,10:90]
#' rpci <- pci/matrix(1,dim(pci))%*%colMeans(pci)
#' m <- mkv(rpci)
#' shorrock(m[[1]])
#'
#' @export

shorrock <- function(x){
	if(ncol(x)!=nrow(x)) stop(x, "is not a square matrix")
	trace <- sum(diag(x))
	n <- nrow(x)
	sh <- (n - trace)/ (n - 1)
	return(sh)
}
