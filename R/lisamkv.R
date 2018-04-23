require(spdep)
require(rgdal)


source("quad.R")


lisamkv<-function(x,W){
  if(is.null(dim(x))==TRUE) stop("You must provide a matrix conteaining n spatial unita and t  time periods")
  t<-dim(x)[2L]
  n<-dim(x)[1L]
  if(t<2L) stop("At least you must provide two time periods")
  LISA <- apply(x,2,quad,W=W)
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

ll <- lisamkv(dat[,9:48],w1rook) 



