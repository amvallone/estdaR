
quad <- function(x,W,geoda=FALSE){
	if(geoda==TRUE){
		qn <- c("hh"=1,"ll"=2,"hl"=4,"lh"=3)
	} else {
		qn <- c("hh"=1,"ll"=3,"hl"=4,"lh"=2)
	}
	zx <- (x - mean(x))#/sd(x)
	lx <- spdep::lag.listw(W,x)
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