#' @import rgdal
#'
#' @name lisamkv
#' @rdname lisamkv
#'
#' @title  Markov for Local Indicators of Spatial Association
#' @description Compute the  Markov transition matrix for Local Indicators of Spatial Association
#'
#' @param x numerical matrix of n spatial unit ans t time periods
#' @param W an objet of listw class.
#' @param ... other argument to \code{quad} function. See \code{\link{quad}} for more infomation.
#' @details For later...
#'
#' @return a list cantaning three object
#' \describe{
#'   \item{move}{a data frame indicating which type of LISA transition occurred}
#'   \item{lisamatrix}{markov LISA transition matrix}
#'   \item{p.lisamatrix}{markov probability LISA transition matrix}
#' }
#'
#' @examples
#' data(us48)
#' data <- as.data.frame(us48)
#' w1queen <- nb2listw(poly2nb(us48))
#' ll <- lisamkv(data[,10:90],w1queen)
#'
#' @export

lisamkv<-function(x,W,...){
  if(is.null(dim(x))==TRUE) stop("You must provide a matrix conteaining n spatial unita and t  time periods")
  t<-dim(x)[2L]
  n<-dim(x)[1L]
  if(t<2L) stop("At least you must provide two time periods")
  LISA <- apply(x,2,quad,W=W,...)
  type <- matrix(c(1:16),nrow=4,ncol=4,byrow=TRUE) #All the possible movements
  move<-matrix(0,nrow=n,ncol=(t-1))
  for (i in 1L:(t-1)){
    for (j in 1:n){
		move[j,i] <- type[ LISA[j,i] , LISA[j,i+1] ]
    }
  }
  f.move <- table(factor(move,levels=c(1:16)))
  lmatrix <- matrix(f.move,nrow=4,ncol=4,byrow=TRUE)
  p.lisa<-lmatrix/rowSums(lmatrix)
  output<-list(move,lmatrix,p.lisa)
  names(output)<-c("move","lisamatrix","p.lisamatrix")
  return(output)
}






