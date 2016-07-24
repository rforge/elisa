# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 1.1
# Licence GPL v3 

.checkrasterMemory <- function(cells,n=1) {
  cells <- ceiling(sqrt(cells))
  canProcessInMemory(raster(nrows=cells, ncols=cells, xmn=0, xmx=cells),n)
}

if (!isGeneric("entrogram")) {
  setGeneric("entrogram", function(x, width, cutoff,...)
    standardGeneric("entrogram"))
}


setMethod('entrogram', signature(x='RasterLayer'), 
          function(x, width, cutoff, categorical, nc, dif, cloud=FALSE, s=NULL,...) {
            re <- res(x)[1]
            if (missing(cutoff)) cutoff<- sqrt((xmin(x)-xmax(x))^2+(ymin(x)-ymax(x))^2) / 3
            if (missing(width)) width <- re
            else if (width < re) width <- re
            if (cutoff < width) stop("cutoff should be greater than width size")
            nlag <- ceiling(cutoff / width)
            
            n <- ncell(x) - cellStats(x,'countNA')
            #---
            if (is.null(s)) {
              if (!.checkrasterMemory(n,nlag)) {
                s <- c()
                for (i in (nlag-1):1) s <- c(s,.checkrasterMemory(n,i))
                s <- which(s)
                if (length(s) > 0) {
                  s <- (nlag - s[1]) / (2*nlag)
                  s <- ceiling(n * s)
                  s <- sampleRandom(x,s,cells=TRUE)[,1]
                } else {
                  s <- 1 / (2 * nlag)
                  s <- ceiling(n * s)
                  while (!.checkrasterMemory(s,1)) s <- ceiling(s / 2)
                  s <- sampleRandom(x,s,cells=TRUE)[,1]
                }
              } else {
                s <- (1:ncell(x))[which(!is.na(x[]))]
              }
            } else {
              if (!is.numeric(s)) stop("s argument should be an integer number or NULL!")
              while (!.checkrasterMemory(s[1],1)) s <- ceiling(s[1] * 0.8)
              if (s > n) s <- n
              s <- sampleRandom(x,s,cells=TRUE)[,1]
            }
            #######---------------
            #----
            if (!missing(nc)) {
              if (missing(categorical)) {
                if (missing(dif)) categorical <- FALSE
                else {
                  categorical <- TRUE
                  cat("input data is considered categorical, and nc is ignored!\n")
                }
              }
            } else {
              if (missing(categorical) && !missing(dif)) categorical <- TRUE
            }
            #----
            if (missing(categorical) || !is.logical(categorical)) {
              # guessing whether the layer is categorical:
              if (.is.categorical(x)) {
                categorical <- TRUE
                cat("the input is considered as a categorical variable...\n")
              } else {
                categorical <- FALSE
                cat("the input is considered as a continuous variable...\n")
              }
            }
            #----
            if (!categorical && missing(nc)) {
              nc <- nclass(x)
            } else if (categorical) {
              classes <- unique(x)
              nc <- length(classes)
            }
            #-----
            
            if (categorical) {
              if (missing(dif)) {
                dif <- rep(1,nc*nc)
                for (i in 1:nc) dif[(i-1)*nc+i] <-0
              } else {
                dif <- .checkDif(dif,classes)
              }
            }
            #-----
            #######---------------------
            out <- new("Entrogram")
            out@width <- width
            out@cutoff <- cutoff
            if (cloud) {
              out@entrogramCloud <- matrix(NA,nrow=length(s),ncol=nlag)
              for (i in 1:nlag) {
                out@entrogramCloud[,i] <- elsa(x,d=i*width,nc=nc,categorical=categorical,dif=dif,cells=s)
              }
              out@entrogram <- data.frame(distance=seq(width,width*nlag,width) - (width/2),E=apply(out@entrogramCloud,2,mean,na.rm=TRUE))
            } else {
              d <- seq(width,width*nlag,width) - (width/2)
              out@entrogram <- data.frame(distance=d,E=rep(NA,length(d)))
              for (i in 1:nlag) {
                out@entrogram [i,2] <- mean(elsa(x,d=i*width,nc=nc,categorical=categorical,dif=dif,cells=s),na.rm=TRUE)
              }
            }
            out
          }
)
