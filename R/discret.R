#' @name discret
#' @rdname discret
#'
#' @title Discretization of a continuous variable
#' @description Transform a continous variable into a discret variable.
#'
#' @param  x numerical vector
#' @param classes a number of a numeric vector of two or more unique cut points giving the number of intervals into which x will be cut
#' @param type 	an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. For more information see the quantile fuction
#' @param ... other argumnt to \code{\link{discret}} function
#'
#' @details for later..
#'
#' @return a vector
#'
#' @examples
#' x <- rnorm(1000)
#' dx <- discret(x,4)
#'
#' @export


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
