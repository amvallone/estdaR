#' @name prais
#' @rdname prais
#'
#' @title Prais conditional mobility measure.
#' @description Compute the Prais conditional mobility measure.
#' @param x Markov probability transition matrix.
#'
#' @details Prais' conditional mobility measure for a class is defined as:
#'     \eqn{ pr_i = 1 -  p_{i,i}}
#'
#' @return a vector
#'
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' pci <- data[,10:90]
#' rpci <- pci/matrix(1,dim(pci))%*%colMeans(pci)
#' m<-mkv(rpci)
#' prais(m[[1]])
#'
#' @export

prais <- function(x){
	pr <- 1-diag(x)
	return(pr)
}



