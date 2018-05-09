#'
#'@name quad
#'@rdname quad
#'
#'@title Quandrant Allocator
#'@description Compute the position of a spatial unit and its neihbors in the MOran Scatter Plot
#'
#'@param x a vector
#'@param W a \code{listw} object.
#'@param y a vector, use it in case that the lag must be compute in another variable.
#'@param geoda logial.  If True use GeoDa quadrant scheme: HH=1, LL=2, LH=3, HL=4. If False use PySAL quadrant Scheme: HH=1, LH=2, LL=3, HL=4
#'
#'@details explain here what is the geoda
#'@return a vector
#'
#'@examples
#'data(Guerry)
#'quad(pc1,w1queen)
#'
#'@export
#'

quad <- function(x,W,y,geoda=TRUE){
  if(missing(y)){y <- x}
	if(geoda==TRUE){
		qn <- c("hh"=1,"ll"=2,"hl"=4,"lh"=3)
	} else {
		qn <- c("hh"=1,"ll"=3,"hl"=4,"lh"=2)
	}
	zx <- (x - mean(x))#/sd(x)
	lx <- spdep::lag.listw(W,y)
	zlx <- (lx - mean(lx))#/sd(lx)
	hx <- zx > 0
	hlx <-  zlx > 0
	hh <- hx * hlx
	lh <- (1 - hx) * hlx
	hl <- hx * (1 - hlx)
	ll <- (1 - hx) * (1 - hlx)
	out <- qn["hh"]*hh + qn["ll"]*ll + qn["hl"]*hl + qn["lh"]*lh
	return(out)
}
