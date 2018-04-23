
# x numerical vector
# classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
# type 	an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. For more information see the quantile fuction 
# internal function used en mkv and sp.mkv functions

discret<-function(x, classes=5,type=7,...){ #subdivide una variable en clases
    if(length(classes)!=1){# class definition ad hoc.
      if(max(classes)==1){ # breaks are probabilities
        aux<-quantile(x,classes,type = type, ...)
        aux<-c(-Inf,aux[2:(length(aux)-1)],Inf)
        output<-as.numeric(cut(x,breaks=aux))
      } else {
        aux<-c(-Inf,classes,Inf)
        output<-as.numeric(cut(x,breaks=aux))
      }
    } else {
      aux <- quantile(x,seq(0,1,1/classes))	
      output<-as.numeric(cut(x,breaks=aux,include.lowest=TRUE))
    }
    return(output)
}