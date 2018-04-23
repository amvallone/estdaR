

join.d <- function(x,W){
	t<-dim(x)[2L]
  	n<-dim(x)[1L]
	x.mean <- colMeans(x)
	x.bar <- x/x.mean
	x.L<-matrix(0,nrow=n,ncol=t)
		for (i in 1:t){
			x.L[,i]<-lag.listw(W,x[,i])
		}
	x.L.mean <- colMeans(x.L)
	x.L.bar <- x.L/(x.L.mean)
	h.x <- apply(x.bar > 1,2,as.numeric) + 1 # 2 are a h and 1 a low
	h.x.l <- apply(x.L.bar > 1,2,as.numeric) + 1
	#################################################
	m.x <- mkv.int(h.x)
	m.x <- m.x/rowSums(m.x)
	m.x.l <- mkv.int(h.x.l)
	m.x.l <- m.x.l/rowSums(m.x.l)
	A <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0),4,4,byrow=TRUE)
	kp <- m.x %x% m.x.l
	aux <- A %*% kp %*% t(A)
	trans <- lisamkv(x,W)$lisamatrix
	t.trans <- rowSums(trans)
	t.h0 <- diag(t.trans) %*% aux
	chi <- m.chi(trans,t.h0)
	output <- list("Chi2"=chi,"Expected"=t.h0)
	return(output)
}


#### Example


us_inc <- readOGR("./exampleData/us_income","us48")
inc_data <- read.csv("./exampleData/us_income/spi_download.csv")
inc_data <- inc_data[which(complete.cases(inc_data)==TRUE),3:43]
us<-inc_data[1,]
inc_data[,1]<-as.character(inc_data[,1])
out=c("United States 3/","Alaska 3/","District of Columbia","Hawaii 3/","New England", "Mideast","Great Lakes","Plains","Southeast","Southwest","Rocky Mountain", "Far West 3/")
p<-c()
for(i in seq_along(out)){ p<-c(p,which(inc_data$AreaName==out[i]))}
inc_data <- inc_data[- p,]
data <- merge(us_inc,inc_data,by.x="STATE_NAME", by.y="AreaName")
rook.nb <- poly2nb(data)
w1rook <- nb2listw(rook.nb,style="W")
dat <- as.data.frame(data)

join.d(dat[,9:48],w1rook) 

