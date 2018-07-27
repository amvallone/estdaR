#' @name homo.test
#' @rdname homo.test
#' 
#' @title Test of time homogeneity
#' 
#' @description Computes the Test of time homogeneity  based on Bickenbach and Bode (2003)
#' 
#' @param x numerical matrix of n spatial unit ans t time periods
#' @param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
#' @param pr number of subperiods in which to divide the entire sample
#' 
#' @details for later.....
#' 
#' @return A list coantaning the Q statistic and the LR statistic.
#' 
#' @references Bickenbach, F. and Bode, E. 2003. Evaluating the Markov Property in Studies of Economic Convergence, International Regional Science Review, vol. 26, no. 3, 363–92
#' 
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' pci <- data[,10:90]
#' rpci <- pci/matrix(1,dim(pci))%*%colMeans(pci)
#' homo.test(rpci,pr=5)
#'
#' @export


homo.test <- function(x,classes=5,pr=3){
  n <- classes
  p <- mkv(x,classes=classes)[[1]]
  A_i <- rowSums(p>0)
  M <- trunc(ncol(x)/pr)
  s <- c(0,c(1L:pr)*M)
  p_ij <- array(0,dim=c(n,n,pr))
  n_ij <- array(0,dim=c(n,n,pr))
  for( i in 1L:(length(s)-1)){
    sub<-x[,(s[i]+1):s[i+1]]
    aux <- mkv(sub,classes)
    p_ij[,,i] <- aux[[1]]
    n_ij[,,i] <- aux[[2]]
  }
  aux <- mkv(x[,(s[length(s)-1]+1):ncol(x)],classes)
  p_ij[,,pr] <- aux[[1]]
  n_ij[,,pr] <- aux[[2]]
  b_i <- rowSums(sapply(seq_len(pr),function(x) rowSums(n_ij[,,x]))>0)
  p<-array(p,dim=c(n,n,pr))
  q1 <- ((p_ij-p)^2)/(p+1*(p==0))
  q2 <- array(unlist(lapply(seq_len(pr),function(x)diag(rowSums(n_ij[,,x]),n,n))),dim=c(n,n,pr))
  q <- sapply(seq_len(pr),function(x) q2[,,x]%*% q1[,,x])
  Q <- sum(q)
  p.s <- (p_ij>0)*(p>0)
  np.s <- 1*(p.s==0)
  l.ratio <- log((p.s*p_ij+np.s)/(p.s*p+np.s))
  LR <- 2*sum(n_ij*l.ratio)
  dof <- sum((A_i-1)*(b_i-1))
  pQ <- 1-pchisq(Q,dof)
  pL.R <- 1-pchisq(LR,dof)
  testQ <- data.frame("Q"=round(Q,2),"p-value"=round(pQ,4),"d.o.f" = round(dof,0))
  testLR <- data.frame("LR"=round(LR,2),"p-value"=round(pL.R,4),"d.o.f"= round(dof,0))
  out <- list("Q"=testQ, "LR"=testLR)
  return(out)
}

