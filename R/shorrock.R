

#' param x Markov probability transition matrix.
shorrock <- function(x){
	if(ncol(x)!=nrow(x)) stop(x, "is not a square matrix")
	trace <- sum(diag(x))
	n <- nrow(x)
	sh <- (n - trace)/ (n - 1)
	return(sh)
}

## example

usinc <- read.csv("./exampleData/us_income/usjoin.csv")
pci <- usinc[,-1]
pci <- pci[,-1]
us_mean<-colMeans(pci)
rpci <- matrix(0,48,81)
for(i in seq_len(ncol(rpci))){rpci[,i] <- pci[,i]/us_mean[i]}
w <- read.gal("./exampleData/us_income/states48.gal",region.id=c(0:47))
w <- nb2listw(w)

m<-mkv(rpci)

shorrock(m[[1]])

