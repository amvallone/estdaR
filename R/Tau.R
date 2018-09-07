#' @import NORMT3
#' @name tau
#' @rdname tau
#' @title Kendall rank correlation coefficient
#' @description Compute the Kendall's tau rank correlation coefficient index
#'
#' @param x rank variable
#' @param y rank variable.
#'
#' @details Kendall’s Tau is a non-parametric measure of relationships between columns of ranked data. The implemtation is based on Rey(2014),
#'
#' @return a vector coantaning \itemize{
#'    \item The Tau statistic value
#'    \item The tau's p value
#'	  \item The number of concordant pairs
#' 	  \item The number of Discordant pairs
#'	  \item The number of extra x pairs. An extra x pair is a pair which \eqn{sgn(x_i - x_j)=0}
#'	  \item The number of extra y pairs. An extra y pair is a pair which \eqn{sgn(y_i - y_j)=0}
#'    }
#' @examples
#' a<-c(12,2,1,12,2)
#' b<-c(1,4,7,1,0)
#' cor(a,b,method="kendall")
#' tau(a,b)
#'
#'
#' @export


tau <- function(x,y){
	if (length(x)!=length(y)) stop("The rank variable must have the same length")
	n <- length(y)
	data <- data.frame("x"=x,"y"=y)
	data <- data[order(data$x,data$y),]
	perm <- as.numeric(rownames(data))
	vals <- y[perm]
	ExtraY <- 0
    ExtraX <- 0
    ACount <- 0
    BCount <- 0
    CCount <- 0
    DCount <- 0
    ECount <- 0
    DCount <- 0
    Concordant <- 0
    Discordant <- 0
    #left child's id
    li <- rep(NA,(n-1))
    #right child's id
    ri <- rep(NA,(n-1))
    # number of left descendants for a node
    ld <- rep(0,n)
    # number of values equal to value i
    nequal <- rep(0,n)

    for (i in 2:n){
    	NumBefore <- 1
    	NumEqual <- 2
    	root <- 1
    	x0 <- x[perm[(i-1)]]
    	y0 <- y[perm[(i-1)]]
    	x1 <- x[perm[i]]
    	y1 <- y[perm[i]]
    	if(x0!=x1){
    		DCount <- 0
    		ECount <- 1
    	} else{
    		if(y0 == y1){
    			ECount <- ECount +1
    		} else {
    			DCount <- DCount + ECount
    			ECount <- 1
    		}
    	}
    	root <- 1
    	inserting <- TRUE
    	while(inserting==TRUE){
    		current <- y[perm[i]]
    		if(current > y[perm[root]]){
    			NumBefore <-  NumBefore + 1 + ld[root] + nequal[root]
    			if(is.na(ri[root])){
    				ri[root] <- i
    				inserting <- FALSE
    			} else {
    				root <- ri[root]
    			}
    		} else if(current < y[perm[root]]){
    			ld[root] <- ld[root]+1
    			if(is.na(li[root])){
    				li[root] <- i
    				inserting <- FALSE
    			} else {
    				root <- li[root]
    			}
    		} else if(current == y[perm[root]]){
    			NumBefore <- NumBefore + ld[root]
    			NumEqual <- NumEqual + nequal[root] +1
    			nequal[root] <- nequal[root] + 1
    			inserting <- FALSE
    		}

    	}
    	ACount <- (NumBefore-1) - DCount # substracting 1 for index different
    	BCount <- (NumEqual-1) - ECount # substracting 1 for index different
    	CCount <- i -(ACount + BCount + DCount + ECount - 1)-1 # substracting 1 for index different
    	ExtraY <- ExtraY + DCount
    	ExtraX <- ExtraX + BCount
    	Concordant <- Concordant + ACount
    	Discordant <-  Discordant + CCount
    }
    cd <- Concordant + Discordant
    num <- Concordant - Discordant
    tau <- num / sqrt((cd + ExtraX) * (cd + ExtraY))
    v <- (4 * n + 10) / (9 * n * (n - 1))
    z <- tau / sqrt(v)
    pval <- NORMT3::erfc(abs(z/1.4142136))
    out <- c("Tau"=tau,"pval"=pval,"Concordant"=Concordant,"Discordant"=Discordant,"ExtraX"=ExtraX,"ExtraY"=ExtraY)
    return(Re(out))
}




