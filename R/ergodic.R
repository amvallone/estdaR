
ergodic <- function(m){
	aux <- eigen(t(m))
	eva <- aux$values
	eve <- aux$vectors
	m.eva <- max(eva)
	i <- which(eva==m.eva)
	return(eve[,i]/sum(eve[,i]))
}

t <- matrix(c(.5,.25,.25,.5,0,.5,.25,.25,.5),3,3,byrow=TRUE)

ergodic(t)

