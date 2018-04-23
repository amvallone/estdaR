

# description Compute the Spatial Markov transition probability matrix
# param x numerical matrix of n spatial unit ans t time periods
# param W an objet of listw class.
# param t period of time to compute the spatial lag.
# param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
# param fixed logical, if it is TRUE the data are pooled over space and time and the quintiles calculated for the pooled data


sp.mkv <- function(m,W,classes=5,fixed=TRUE){
	lag<-matrix(0,nrow(m),ncol(m))
	for(j in 1:ncol(m)){
		lag[,j]<-lag.listw(W,m[,j])
	}
	if(fixed==TRUE){
		auxl <- unlist(c(lag))
		auxx <- unlist(c(m))
		auxl <- discret(auxl,classes=classes)
		auxx <- discret(auxx,classes=classes)
		lx <- matrix(auxl,dim(lag))
		x <- matrix(auxx,dim(m))
	} else {
		lx <- apply(lag,2,discret,classes=classes) 
		x <- apply(m,2,discret,classes=classes)
	}
	class<-as.numeric(unique(as.factor(lx)))
    class<-sort(class)
    mm <- array(0,dim = c(length(class),length(class),length(class)),dimnames=list(class,class,paste("Lag",class)))
	for (j in 1:ncol(x)-1){
		for(i in 1:nrow(x)){
        		mm[x[i,j], x[i,j+1],lx[i,j]] <- mm[x[i,j], x[i,j+1],lx[i,j]] + 1
      		}
    }
    mm.p <- array(0,dim = c(length(class),length(class),length(class)),dimnames=list(class,class,paste("Lag",class)))
    for(k in seq_along(class)){
    		mm.p[,,k] <- mm[,,k]/rowSums(mm[,,k])
    }
  output <- list("Probabilities" = mm.p, "Transitions" = mm)
  return(output)
 }


#example
usinc <- read.csv("./exampleData/us_income/usjoin.csv")
pci <- usinc[,-1]
pci <- pci[,-1]
us_mean<-colMeans(pci)
rpci <- matrix(0,48,81)
for(i in seq_len(ncol(rpci))){rpci[,i] <- pci[,i]/us_mean[i]}
w <- read.gal("./exampleData/us_income/states48.gal",region.id=c(0:47))
w <- nb2listw(w)




# Spatial Markov
sm <- sp.mkv(rpci,w)
