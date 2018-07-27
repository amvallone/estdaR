#' @name sig.lisamkv
#' @rdname sig.lisamkv
#'
#' @title  Markov for Local Indicators of Spatial Association including a significant state
#' @description Compute the  Markov transition matrix for Local Indicators of Spatial Association based in the significance of the local moran indicator
#'
#' @param x numerical matrix of n spatial unit ans t time periods
#' @param W an objet of listw class.
#' @param nsim  number of random spatial permutations for calculation of pseudo p-values, by default is set in 999
#' @param ... other argument to \code{quad} function. See \code{\link{quad}} for more infomation.
#' @details first the the Local Moran indicator is compute useing the \code{nsim} number of permutation. The non significant indication is considere as an state of the markow chain.
#'
#' @return a list cantaning three object
#' \describe{
#'   \item{move}{a data frame indicating which type of LISA transition occurred}
#'   \item{lisamatrix}{markov LISA transition matrix}
#'   \item{p.lisamatrix}{markov probability LISA transition matrix}
#' }
#'
#' @examples
#' \dontrun{
#' data(us48)
#' data <- as.data.frame(us48)
#' w1queen <- nb2listw(poly2nb(us48))
#' ll <- sig.lisamkv(data[,10:90],w1queen,999,geoda=FALSE)
#'}
#' @export

sig.lisamkv <- function(x,W,nsim=999,...){
  if(is.null(dim(x))==TRUE) stop("You must provide a matrix conteaining n spatial unita and t  time periods")
  t<-dim(x)[2L]
  n<-dim(x)[1L]
  if(t<2L) stop("At least you must provide two time periods")
  l1 <- lapply(1:t,function(i) moran(x[,i],W,nsim=nsim,...))
  states <- matrix(0,n,t)
  for (i in 1:t){
    states[,i] <- l1[[i]][,6L]
  }
  states <- ifelse(states==0,5,states)
  type <- matrix(c(1:25),nrow=5,ncol=5,byrow=TRUE) #All the possible movements
  move<-matrix(0,nrow=n,ncol=(t-1))
  for (i in 1L:(t-1)){
    for (j in 1:n){
      move[j,i] <- type[ states[j,i] , states[j,i+1] ]
    }
  }
  f.move <- table(factor(move,levels=c(1:25)))
  lmatrix <- matrix(f.move,nrow=5,ncol=5,byrow=TRUE)
  p.lisa<-lmatrix/rowSums(lmatrix)
  output<-list(move,lmatrix,p.lisa)
  names(output)<-c("move","lisamatrix","p.lisamatrix")
  return(output)
}