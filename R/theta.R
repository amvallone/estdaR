

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

Reg <- c(rep(1,4),rep(2,4),rep(3,4),rep(4,4))
x0 <- c(1:16)
x1 <- c(13:16,5:12,1:4)
x11 <- c(4:1,5:16)

theta.int(x0,x1,Reg)


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



mexico <- read.csv("./exampleData/mexico.csv")

theta(mexico[,1:7],mexico[,8],999)