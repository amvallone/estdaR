#' @import gridExtra RColorBrewer spdep ggplot2 stats graphics
#'
#'
#' @name d.LISA
#' @rdname d.LISA
#'
#' @title Directional LISA
#' @description Compute the origin standardized movement of a spatial unit and its neighbors and the pseudo p-value of the directional co-movement.
#'
#' @param x0 a vector conatining the variale at initial period of analyis
#' @param x1 a vector conatining the variale at final period of analyis
#' @param W an objecto of class \code{listw}
#' @param Regime a numerci vector contaning the regime to the sptail unit belongs
#' @param k a scalar in c(4,8) indicating number of circular sectors in rose diagram, by default it is set as 8
#' @param mean.rel a logical, Is the data mean relative?. By default it is FALSE
#' @param nsim  number of random spatial permutations for calculation of pseudo p-values, the default value is NULL.
#'
#' @details For later...
#'
#' @references Rey, S. J., Murray, A. T., & Anselin, L. (2011). Visualizing regional income distribution dynamics. Letters in Spatial and Resource Sciences, 4(1), 81–90.
#'
#' @return a list contaning \describe{
#'     \item{"Lisa"}{A vector scatterplot showing the Directinal co-movent of a spatial unit and its neighbors}
#'     \item{"rose"}{A rose diagram of the LISA}
#'     \item{"p.rose"}{pseudo p-values of the direction in the rose diagram. Only available when nsimm is not NULL}
#'     \item{"data"}{a data frame contening the data used for the lisa scatterplot}
#'     \item{"counts"}{a data frame conating the information used for the rose diagram}
#'     \item{"p.value"}{a data frame contaning the information used in the pseudo p-value graph. ONly available when nsim is not NULL}
#'  }
#'
#' @examples
#' data(us48)
#' t0<-data$X1969../us.total$X1969..
#' t1 <-data$X2008../us.total$X2008..
#' Regime <-data$SUB_REGION
#' ok <- d.LISA(t0,t1,w1rook,Regime,k=8,nsim=999)
#'
#' @export



d.LISA <- function(x0,x1,W,Regime=NULL,k=8,mean.rel=FALSE,nsim=NULL){
	skip.leg <- ifelse(is.null(Regime)==TRUE,1,0)
	if (is.null(Regime)==TRUE) { Regime <- rep(1,length(x0)) }
	x0.lag <- lag.listw(W,x0)
	x1.lag <- lag.listw(W,x1)
	if (mean.rel==TRUE){
		x0 <- x0 - mean(x0)
		x1 <- x1 - mean(x1)
		x0.lag <- x0.lag - mean(x0.lag)
		x1.lag <- x1.lag - mean(x1.lag)
	}
	Regime <- factor(Regime)
	f.point <- as.data.frame(cbind(x1-x0,x1.lag-x0.lag))
	vtop<-cbind(rep(0,length(x0)),rep(0,length(x0)),f.point,Regime)
	colnames(vtop)<-c("x0","y0","x1","y1")
	lisa <-ggplot(vtop)+
		geom_hline(yintercept=0,linetype=1)+
		geom_vline(xintercept=0,linetype=1)+
		geom_segment(aes(x=x0, xend=x1, y=y0, yend=y1,color=Regime),size=1.1,arrow=arrow(length = unit(0.2, "cm")))+
		theme(
		legend.position="bottom",
		panel.border = element_rect(linetype = "solid", fill = NA),
		panel.background = element_rect(fill = NA),
		panel.grid.major = element_line()
			)+xlab("X")+ylab("WX")+scale_color_brewer(palette="Set1",direction=-1)

#_________________ROSE_____________________________________
	rotar <-ifelse(k==8,pi/4,0)
	step <- (2*pi)/k
	breaks <- seq(0,2*pi,step)
	lmt <- (breaks * 180)/pi # rad to deg
	symb<-c()
	for (i in seq_len(length(lmt)-1)){
		symb <- c(symb, paste(lmt[i],lmt[i+1],sep="-"))
	}
	z<-atan2(f.point[,2],f.point[,1]) #get angles
	z <- ifelse(z>0,1*z,(2*pi)+z) #avoid negatives rad
	bins = rep(1, length(z))
    i = 1L
    for(b in breaks){
        bins[z > b] = i
        i = i + 1L
    }
	angle <- factor(bins, levels = seq_along(symb), labels = symb)
	len<-c()
	for (i in seq_len(length(x0))){
		len<-c(len,spDistsN1(matrix(c(vtop[i,1L],vtop[i,2L]),nrow=1),matrix(c(vtop[i,3L],vtop[i,4L]),nrow=1)))
	}
	real<-hist(z,breaks=breaks,plot=FALSE)$counts
	d.plot <- data.frame(Regime=factor(Regime),angle=as.character(angle),length=len)
	d.plot <- within(d.plot,angle <- factor(angle,levels=c(symb[1],symb[length(lmt):2]))) # inverse order
	rose<-ggplot(d.plot,aes(angle,fill=Regime),xlab=" ",ylab=" ") +geom_bar(width=1,colour="black",size=0.1) +scale_x_discrete(name="")+scale_y_continuous(name="",breaks=seq(0,max(table(angle)),5),labels=seq(0,max(table(angle)),5))+coord_polar(start=rotar)+theme(legend.position="bottom",panel.grid.major = element_line( color="gray50",linetype="solid"),panel.background = element_blank())+scale_fill_brewer(palette="Set1",direction=-1)

#__________P VALUES_______________________

	if(!is.null(nsim)){
	id <- seq_along(x0)
	capture <- matrix(0,nsim,k)
	for (m in seq_len(nsim)){
		rid <- sample(id)
		rx0<-x0[rid]
		rx1<-x1[rid]
		rx0.lag <- lag.listw(W,rx0)
		rx1.lag <- lag.listw(W,rx1)
		rf.point <- as.data.frame(cbind(rx1-rx0,rx1.lag-rx0.lag))
		rz<-atan2(rf.point[,2],rf.point[,1])
		rz <- ifelse(rz>0,1*rz,(2*pi)+rz)
		rf<-hist(rz,breaks=breaks,plot=FALSE)$counts
		capture[m,]<-rf
	}
	colnames(capture) <-symb
	larger <- c()
		for (t in seq_along(real)){
			larger <- c(larger, sum(capture[,t]>=real[t]))
		}
	p.l <- nsim - larger
	p <- ifelse(p.l<larger, p.l,larger)
	p.value <-(p+1)/(nsim+1)
	p.plot<-rep(1,length(p.value))
	p.plot[p.value<=0.001]<-0.001
	p.plot[p.value>0.001 & p.value<=0.01]<-0.01
	p.plot[p.value>0.01 & p.value<=0.05]<-0.05
	#p.value[p.value>0.05 & p.value<=0.1]<-0.1
	p.plot<-factor(p.plot)
	d.pv <- data.frame(type=symb,p.value=p.plot)
	expec <- as.vector(sapply(as.data.frame(capture),mean))
	d.pv <- within(d.pv,type <- factor(type,levels=c(symb[1],symb[length(lmt):2])))
	val <- c("0.001" = "darkblue", "0.01" = "steelblue4","0.05"="skyblue2", "1"="azure2")
	lab <- c("0.001", "0.01", "0.05", "1")
	are <- which(lab %in% p.plot)
	pv <- ggplot(d.pv,aes(type,fill=p.value))+geom_bar(width=1)+coord_polar(start=rotar)+theme(panel.grid.major = element_line( color="gray50",linetype="solid"),panel.background = element_blank())+ylab("")+xlab("")+scale_y_continuous(breaks=NULL)+ scale_fill_manual(values = val,labels = lab ,limits = lab)
		}

#________ OUTPUT_________________________

	if(skip.leg==1){
		lisa <- lisa + theme(legend.position="none")
		rose <- rose + theme(legend.position="none")
	} else {
		lisa <- lisa #+ scale_linetype_discrete(name="Regime")
		rose <- rose + theme(legend.position="right")
	}
	if(!is.null(nsim)){
		d.real<- data.frame(type=symb,Counts=real,Expected=expec)
		gridExtra::grid.arrange(lisa,arrangeGrob(rose, pv, nrow = 2), ncol = 2)
		output<-list(lisamap=lisa,rose=rose,p.rose=pv,data=d.plot,counts=d.real,p.value=d.pv)
	}else{
		d.real<- data.frame(type=symb,Counts=real)
		gridExtra::grid.arrange(lisa,rose,nrow=1,ncol=2)
		output<-list(lisamap=lisa,rose=rose,data=d.plot,counts=d.real)
	}
	invisible(output)
}



