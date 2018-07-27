

d.lisamkv<-function(x0,x1,W,mean.rel=FALSE,arrow=TRUE,...){
  q.x0 <- quad(x0,W,...)
  q.x1 <- quad(x1,W,...)
  type <- matrix(c(1:16),nrow=4,ncol=4,byrow=TRUE) #All the possible movements
  move <- diag(type[q.x0,q.x1])
  x0.lag <- lag.listw(W,x0)
  x1.lag <- lag.listw(W,x1)
  if (mean.rel==TRUE){
    x0 <- x0 - mean(x0)
    x1 <- x1 - mean(x1)
    x0.lag <- x0.lag - mean(x0.lag)
    x1.lag <- x1.lag - mean(x1.lag)
  }
  f.point <- as.data.frame(cbind(x1-x0,x1.lag-x0.lag))
  y0 <- NULL # unuseful in the code, but avoids a CRAN note.
  y1 <- NULL
  vtop<-cbind(rep(0,length(x0)),rep(0,length(x0)),f.point,"move"=factor(move))
  colnames(vtop)<-c("x0","y0","x1","y1")
  if(arrow==FALSE){
    lisa <-ggplot(vtop)+
      geom_hline(yintercept=0,linetype=1)+
      geom_vline(xintercept=0,linetype=1)+
      geom_point(aes(x=x1, y=y1,color=as.factor(move)),size=2)+
      theme(
        legend.position="bottom",
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line()
      )+xlab("X")+ylab("WX")
  } else {
    lisa <-ggplot(vtop)+
      geom_hline(yintercept=0,linetype=1)+
      geom_vline(xintercept=0,linetype=1)+
      geom_segment(aes(x=x0, xend=x1, y=y0, yend=y1,color=as.factor(move)),size=1.1,arrow=arrow(length = unit(0.2, "cm")))+
      theme(
        legend.position="bottom",
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line()
      )+xlab("X")+ylab("WX")
  }
  lisa
}