#' @name ergodic
#' @rdname ergodic
#'
#' @title Ergodic distribution of a Marov's probability transition matrix
#' @description Calulates the ergodic distribution of a Marov's probability transition matrix
#'
#' @param m Markov probability transition matrix
#'
#' @details for later...
#'
#' @return a vector coantaning the ergodic distribution.
#'
#' @examples
#' t <- matrix(c(.5,.25,.25,.5,0,.5,.25,.25,.5),3,3,byrow=TRUE)
#' ergodic(t)
#'
#' @export
#'

ergodic <- function(m){
	aux <- eigen(t(m))
	eva <- aux$values
	eve <- aux$vectors
	m.eva <- max(eva)
	i <- which(eva==m.eva)
	return(eve[,i]/sum(eve[,i]))
}



