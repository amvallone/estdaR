require(spdep)


# description Compute thre Markov transition probability matrix
# param m numerical matrix of n spatial unit ans t time periods
# param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
# param fixed logical, if it is TRUE the data are pooled over space and time and the quintiles calculated for the pooled data
# param type an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. For more information see the quantile fuction 




mkv<-function(m,classes=5,fixed=FALSE,type=7,...){ #Compute Markox transition matrix
	#Argument checks
  if(is.null(dim(m))==TRUE) stop("You must provide a matrix conteaining n spatial unita and t periods of time")
  t<-dim(m)[2]
  n<-dim(m)[1]
  if(t<2) stop("At least you must provide two period of time for this analysis")
  # pool or not cuts.
  if (fixed==TRUE){
    if(length(classes)!=1){
      if(max(classes)==1){
        cuts<-quantile(as.vector(m),classes,type = type, ...)
      }
    } else {
      cuts<-quantile(as.vector(m) ,seq(0,1,1/5),type = type, ...)
    }
    x<-apply(m,2,discret,classes=cuts)
  } else{
    x<-apply(m,2,discret,classes=classes)
  }
  mm<-mkv.int(x)
  #rownames(mm)<-clases
  #colnames(mm)<-clases
  mm.p <- mm/rowSums(mm)
  output <- list("Probabilities" = mm.p, "Transitions" = mm)
  return(output)
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
