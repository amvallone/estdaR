#' @name geary
#' @rdname guery
#'
#' @title Mutivariable and univariable Geary's C statistic
#' @description Compute the univariable Anselin(1995) and multivarible Anselin(2017) Geary's C statistic
#'
#' @param x a ventor, matrix or data frame containing the variables and spatial units
#' @param W a \code{listw} object
#' @param nsim number of random permutations used to the compute the pseudo p value. By default it is NULL
#' @param type a character indicating the Geary's C to be compute. If it is set as "uni" the standart Geary's C statistic is calulated. Set type as "multi" to cumpute the multivariable Geary's C.
#' @param nbcom number of comparisons use in the Bonferroni multiple comparisons correction. By default it is set aas the number of spatial unit to use.
#' @param ... other argument to \code{quad} function. See \code{\link{quad}} for more information.
#' @details later...
#' @return a data frame (more explanation in Details later)
#' @references
#' Anselin, L. (1995). Local Indicators of Spatial Association—LISA. Geographical Analysis, 27(2), 93–115. \url{https://doi.org/10.1111/j.1538-4632.1995.tb00338.x}
#' Anselin, L  (2017). A Local Indicator of Multivariate Spatial Association: Extending Geary’s c. Center for Spatial Data Science working papers. University of Chicago.
#'
#' @examples
#' data(Guerry)
#' w1queen<-nb2listw(poly2nb(Guerry))
#' pcc <- cbind(Guerry$pc1,Guerry$pc2)
#' e1 <- geary(Guerry$pc1,w1queen,nsim=999)
#' e2 <- geary(pcc,w1queen,nsim=999)
#' e3 <- geary(pcc,w1queen,type="multi",nsim=999)
#' @export


geary <- function(x,W,nsim=NULL,type="uni",nbcom=length(W$neighbours),...){
  #argum check


  # dim check
  if(is.null(dim(x))==TRUE){
    out <- local.c(x,W=W,nsim=nsim,nbcom=nbcom,...)
  } else {
    switch(type,
           uni={
             y <- lapply(seq_len(ncol(x)), function(i) x[,i])
             names(y)<-colnames(x)
             out <- lapply(y,local.c,W=W,nsim=nsim,nbcom=nbcom,...)
           },
           multi={
             out <- m.local.c(x,W=W,nsim=nsim,nbcom=nbcom)
           }
    )
  }
  return(out)
}


# Univariate local Geary

local.c<-function(x,W,nbcom=length(W$neighbours),nsim=NULL,...){
	 l.c<-rep(0,length(x))
	 p.val<-rep(0,length(x))
	 p.l.c<-rep(0,length(x))
	 m2<-var(x)
	 for (i in seq_along(x)){
	 	neg<-W$neighbours[[i]]
	 	w<-W$weights[[i]]
	 	d<-x[i]-x[neg]
	 	c1<-w*(d^2)
	 	c1<-(1/m2)*sum(c1)
	 	l.c[i]<-c1
	 	if(!is.null(nsim)){
	 		p<-rep(0,nsim)
	 		for (j in seq_len(nsim)){
	 			if (i==1) { aux <- c(x[i],sample(x[2:length(x)])) }
	 			 else if (i==length(x)){
	 				aux <- c(sample(x[1:(length(x)-1)]),x[i])
	 			} else{
	 				aux1 <-sample(c(x[1:(i-1)],x[(i+1):length(x)]))
	 				aux<- c(aux1[1:(i-1)],x[i],aux1[i:length(aux1)])
	 				}
	 			d<-x[i]-aux[neg]
	 			c2<-w*(d^2)
	 			c2<-(1/m2)*sum(c2)
	 			p[j]<-c2
	 			}
	 		#plot(density(p))
	 		#abline(v=c1)
	 		above <- p>=c1
	 		larger <- sum(above)
	 		lower <- (nsim - larger)
	 		final <- pmin(larger,lower)
	 		p.val[i]<- (final+1)/(nsim+1)
	 		}
	 	}
	if(!is.null(nsim)){
		p.l.c <-  p.adjust(p.val,method="bonferroni",n=nbcom)
	 	p.plot <- rep(1,length(p.val))
  	p.plot[p.val>=0.01 & p.val<=0.05] <- 0.05
  	p.plot[p.val>0.001 & p.val<=0.01] <- 0.01
		p.plot[p.val<=0.001] <- 0.001
		quad.i <- quad(x,W,...)
		cluster <- quad.i*as.integer(p.plot!=1)
	 	out<-cbind("local.c"=l.c,"p.values"=p.val,"adj. pvalue"=p.l.c,"p Map"=p.plot,"Cluster Map"=cluster)
	} else{
		out<-cbind("local.c"=l.c)
		}
	return(out)
}

# Multivariate Local Geary

m.local.c <- function(X,W,nsim=NULL,nbcom=length(W$neighbours)){
	if(is.data.frame(X)) X<-as.matrix(X)
	if(dim(X)[2L]<2L) stop("You must provide at least 2 variables")
	lc.x <-apply(X,2,local.c,W=W)
	m.lc <- apply(lc.x,1,mean)
	if(!is.null(nsim)){
		p<-matrix(0,nrow(X),nsim)
		larger <- rep(0,nrow(X))
		id <- seq_len(nrow(X))
		for (t in id){
			for (j in seq_len(nsim)){
				if (t==1) {
					orden <- c(id[t],sample(id[2L:length(id)]))
				} else if (t==length(id)){
	 				orden <- c(sample(id[1L:(length(id)-1L)]),id[t])
	 			} else{
	 				aux1 <-sample(c(id[1L:(t-1L)],id[(t+1L):length(id)]))
	 				orden<- c(aux1[1L:(t-1L)],id[t],aux1[t:length(aux1)])
	 			}
	 		X1 <- X[orden,]
	 		rml.c <- c()
	 		for (k in seq_len(ncol(X))){
	 			neg<-W$neighbours[[t]]
	 			w<-W$weights[[t]]
	 			d<-X1[t,k]-X1[neg,k]
	 			c1<-w*(d^2)
	 			c1<-(1/var(X1[,k]))*sum(c1)
	 			rml.c <- c(rml.c,c1)
	 			}
	 		p[t,j]<-mean(rml.c)
			}
		larger[t] <- sum(p[t,]>=m.lc[t])
		lower <- (nsim - larger[t])
		larger[t] <- pmin(larger[t],lower)
		}
		p.value <- (larger+1)/(nsim+1)
		p.adj <- p.adjust(p.value,method="bonferroni",n=nbcom)
		p.plot <- rep(1,length(p.value))
    p.plot[p.value>=0.01 & p.value<=0.05]<-0.05
   	p.plot[p.value>0.001 & p.value<=0.01]<-0.01
    p.plot[p.value<=0.001]<-0.001
    cluster <- as.integer(p.plot!=1)
		out<-cbind("multi local c"=m.lc,"p.value"=p.value,"adj pvalue"=p.adj,"p Map"=p.plot,"Cluster Map"=cluster)
	} else {
		out<-cbind("multi local c"=m.lc)
	}
	return(out)
}





