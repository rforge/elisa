# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

if (!isGeneric("plot")) {
  setGeneric("plot", function(x,y,...)
    standardGeneric("plot"))
}	


setMethod("plot", signature(x='Entrogram'), 
          function(x,xlim,ylim,xlab,ylab,pch,col,main,cloud=FALSE,box=FALSE,...) {
            if ((cloud || box) & is.null(x@entrogramCloud)) stop("cloud=TRUE, or box=TRUE is not working as the cloud data are not stored in the input Entrogram object, you should use cloud=TRUE in the entrogram function to keep the data!")
            nlag <- ceiling(x@cutoff / x@width)
            if (missing(xlim)) xlim <- c(0,x@cutoff)
            if (missing(ylim)) {
              if (cloud || box) ylim <- c(0,quantile(x@entrogramCloud,prob=0.99,na.rm=TRUE))
              else ylim <- c(0,max(x@entrogram$E,na.rm=TRUE))
            }
            if (missing(xlab)) xlab <- "Distance"
            if (missing(ylab)) ylab <- "Mean ELSA"
            if (missing(pch)) pch <- 16
            if (missing(col)) {
              if (box & !cloud) col <- 0
              else col <- 'blue'
            }
            if (missing(main)) {
              if (cloud) main <- "Entrogram Cloud"
              else if (box) main <- "Box plot of Entrogram Cloud"
              else main <- "Entrogram"
            }
            if (cloud) {
              plot(x@entrogram$distance,x@entrogramCloud[1,],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,pch=pch,col=col,...)
              for (i in 2:nlag) points(x@entrogram$distance,x@entrogramCloud[i,],col=col,pch=pch,...)
            } else if (box) boxplot(x@entrogramCloud,names=x@entrogram$distance,xlab=xlab,ylab=ylab,ylim=ylim,col=col,main=main,...)
            else plot(x@entrogram$distance,x@entrogram$E,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,pch=pch,col=col,...)
          }
)

