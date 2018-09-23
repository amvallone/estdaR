#' @name moran
#' @rdname moran
#'
#' @title Univariable and Bivariable Local Moran's I
#' @description Compute the Univariable (Anselin,1995) and Bivariable (cite??) Local Moran's I
#'
#'
#' @param x a vector, matrix or data frame containing the variables and spatial units
#' @param W a \code{listw} object
#' @param nsim number of random permutation used to the compute the pseudo p value. By default it is NULL
#' @param type a character indicating the local Moran's I to be compute. If it is set as "uni" the standard local Moran's I statistic is calulated. Set type as "multi" to cumpute the bivariable local Moran's I.
#' @param nbcom number of comparisons use in the Bonferroni multiple comparisons correction. By default it is set as the number of spatial unit to use.
#' @param ... other argument to \code{quad} function. See \code{\link{quad}} for more infomation.
#' @details later...
#' @return a data frame or a list of data frames (more explanation in Details later)
#' @references
#' Anselin, L. (1995). Local Indicators of Spatial Association—LISA. Geographical Analysis, 27(2), 93–115. \url{https://doi.org/10.1111/j.1538-4632.1995.tb00338.x}
#'
#' @examples
#' data(Guerry)
#' w1queen<-nb2listw(poly2nb(Guerry))
#' pcc <- cbind(Guerry$pc1,Guerry$pc2)
#'\dontrun{ e4 <- moran(Guerry$pc1,w1queen,nsim=999)
#' e6 <- moran(pcc,w1queen,type="multi",nsim=999,geoda=FALSE)}
#'
#' @export


moran <- function(x,W,nsim=NULL,type="uni",nbcom=length(W$neighbours),...){
  #argum check

  # dim check
  if(is.null(dim(x))==TRUE){
    out <- local.m(x,W=W,nsim=nsim,nbcom=nbcom,...)
  } else {
    y <- lapply(seq_len(ncol(x)), function(i) x[,i])
    names(y)<-colnames(x)
    switch(type,
           uni={
             y <- lapply(seq_len(ncol(x)), function(i) x[,i])
             names(y)<-colnames(x)
             out <- lapply(y,local.m,W=W,nsim=nsim,nbcom=nbcom,...)
           },
           multi={
             pairs <- expand.grid(c(1:ncol(x)),c(ncol(x):1))
             pairs <- pairs[which(pairs[,1]!=pairs[,2]),] #avoid the univariable computation.
             out <- list()
             name<-c()
             for( i in seq_len(nrow(pairs))){
               name <- c(name,paste(colnames(x)[pairs[i,1]],colnames(x)[pairs[i,2]],sep="-"))
               m.I<-b.local.m(x[,pairs[i,1]],x[,pairs[i,2]],W=W,nsim=nsim,nbcom=nbcom,...)
               out[[i]] <- m.I
             }
             names(out) <- name
           }
    )
  }
  return(out)
}


# Univariate local Moran

local.m<-function(x,W,nsim=NULL,nbcom=length(W$neighbours),...){
  l.m <- rep(0,length(x))
  p.val <- rep(0,length(x))
  p.l.m <- rep(0,length(x))
  meanx <- mean(x)
  wx <- lag(W,x)
  m.wx <- wx - mean(wx)
  m.x <- x-meanx
  stdx <- sd(x)
  quad <- quad(x,W,...)
  for (i in seq_along(x)){
    neg<-W$neighbours[[i]]
    w<-W$weights[[i]]
    m1 <- w*((x[neg]-meanx)/stdx)
    m1 <- ((x[i]-meanx)/stdx)*sum(m1)
    l.m[i]<-m1
    if(!is.null(nsim)){
    p<-rep(0,nsim)
    for (j in seq_len(nsim)){
      if (i==1) { aux <- c(x[i],sample(x[2:length(x)])) }
      	if (i==length(x)){
        		aux <- c(sample(x[1:(length(x)-1)]),x[i])
      	} else{
        		aux1 <-sample(c(x[1:(i-1)],x[(i+1):length(x)]))
        		aux<- c(aux1[1:(i-1)],x[i],aux1[i:length(aux1)])
      		}
      	m2 <- w*((aux[neg]-meanx)/stdx)
      	m2 <- ((x[i]-meanx)/stdx)*sum(m2)
      	p[j]<-m2
    		}
    	#	plot(density(p))
    	#	abline(v=m1)
		above <- p>=m1
		larger <- sum(above)
		lower <- (nsim - larger)
		final <- pmin(larger,lower)
		p.val[i]<- (final+1)/(nsim+1)
  	}
  }
  if(!is.null(nsim)){
  	p.l.m <-  p.adjust(p.val,method="bonferroni",n=nbcom)
  	p.plot <- rep(1,length(p.val))
  	p.plot[p.val>=0.01 & p.val<=0.05]<-0.05
 	  p.plot[p.val>0.001 & p.val<=0.01]<-0.01
  	p.plot[p.val<=0.001]<-0.001
  	cluster <- quad*as.integer(p.plot!=1)
  	out<-cbind("local m"=l.m,"p.value"=p.val,"adj. p.value"=p.l.m, "Moran quad"=quad,"p Map"=p.plot,"Cluster Map"=cluster)
  	} else {
  		out <- cbind("local.m"=l.m,"Moran quad"=quad)
  	}
  return(out)
}

#Multivariable local Moran

b.local.m<-function(x,y,W,nsim=NULL,nbcom=length(W$neighbours),...){
  l.m <- rep(0,length(x))
  p.val <- rep(0,length(x))
  p.l.m <- rep(0,length(x))
  wy <-lag(W,y)
  m.x <- x - mean(x)
  m.wy <- wy - mean(wy)
  zx <- m.x/sd(x)
  zy <- (y-mean(y))/sd(y) #(mean(y) - y)/sd(y) in geoda
  quad <- quad(x,y=y,W=W,...)
  for (i in seq_along(x)){
    neg<-W$neighbours[[i]]
    w<-W$weights[[i]]
    m1 <- zx[i]*sum(w*zy[neg])
    l.m[i] <- m1
    if(!is.null(nsim)){
    		p<-rep(0,nsim)
    	for (j in seq_len(nsim)){
      		if (i==1) { aux <- c(y[i],sample(y[2:length(y)])) }
      		if (i==length(x)){
        		aux <- c(sample(y[1:(length(x)-1)]),y[i])
      		} else{
        		aux1 <-sample(c(y[1:(i-1)],y[(i+1):length(y)]))
        		aux<- c(aux1[1:(i-1)],y[i],aux1[i:length(aux1)])
      		}
      		zaux <-(aux - mean(aux))/sd(aux) # (mean(aux) - aux)/sd(aux)
      		m2 <- zx[i]*sum(w*zaux[neg])
      		p[j]<-m2
    	}
    #	plot(density(p))
    #	abline(v=m1)
  		above <- p>=m1
	  	larger <- sum(above)
  		lower <- (nsim - larger)
  		final <- pmin(larger,lower)
  		p.val[i]<- (final+1)/(nsim+1)
  		}
    }
	if(!is.null(nsim)){
		p.l.m <-  p.adjust(p.val,method="bonferroni",n=nbcom)
  	p.plot <- rep(1,length(p.val))
  	p.plot[p.val>=0.01 & p.val<=0.05]<-0.05
  	p.plot[p.val>0.001 & p.val<=0.01]<-0.01
  	p.plot[p.val<=0.001]<-0.001
  	cluster <- quad*as.integer(p.plot!=1)
 		out<-cbind("local m"=l.m,"p.value"=p.val,"adj. p.value"=p.l.m, "Moran quad"=quad,"p Map"=p.plot,"Cluster Map"=cluster)
	} else {
		out<-cbind("local m"=l.m, "Moran quad"=quad)
	}
  return(out)
}



