# compute the markov transition matrix
# x is a matrix of discret variables
# internal function used in mkv and join.d

mkv.int <- function(x){
	clases<-as.numeric(unique(as.factor(x)))
  clases<-sort(clases)
	mm<-matrix(0,nrow=length(clases),ncol=length(clases))
    for (i in 1:dim(x)[1]){
      for(j in 1:(dim(x)[2]-1)){
        mm[x[i,j], x[i,j + 1]] <- mm[x[i,j], x[i,j+1]] + 1
      }
    }
    mm
}
