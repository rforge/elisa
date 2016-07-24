# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

if (!isGeneric("elsa.test")) {
  setGeneric("elsa.test", function(x, d, n, nc, categorical, dif,...)
    standardGeneric("elsa.test"))
}


setMethod('elsa.test', signature(x='RasterLayer'), 
          function(x, d, n=99, nc, categorical, dif,cells,filename,...) {
            
            if (missing(filename)) filename <- ''
            
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
            if (!categorical) {
              if (missing(nc)) nc <- nclass(x)
              classes <- 1:nc
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
            
            nNA <- which(!is.na(x[]))
            null.model <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
            null.model <- null.model[nNA]
            if (missing(cells)) {
              e1 <- elsa(x,d=d,nc=nc,categorical=categorical,dif=dif)
              o2 <- x
              #o2 <- calc(o2,function(x) { x[!is.na(x)] <- 0; x})
              o2[nNA] <- 0
              o1 <- raster(x)
              for (i in 1:n) {
                #o1 <- calc(x,function(x) { x[!is.na(x)] <- sample(null.model,length(x[!is.na(x)]),replace=TRUE); x})
                o1[nNA] <- sample(null.model,length(null.model),replace=TRUE)
                e2 <- elsa(o1,d1=d1,d2=d2,nc=nc,categorical=categorical,dif=dif)
                ee <- e1 - e2
                ee <- calc(ee,function(x) {x[x > 0] <- 1; x[x <= 0] = 0; x})
                o2 <- o2 + ee
              }
              rm(e1,e2,ee,o1)
              o2 <- (o2+1) / (n+1)
              
              filename <- trim(filename)
              if (filename != '') writeRaster(o2,filename,...)
            } else {
              e1 <- elsa(x,d=d,nc=nc,categorical=categorical,dif=dif,cells=cells)
              o2 <- rep(0,length(cells))
              for (i in 1:n) {
                o1[nNA] <- sample(null.model,length(null.model),replace=TRUE)
                #o1 <- calc(x,function(x) { x[!is.na(x)] <- sample(classes,length(x[!is.na(x)]),replace=TRUE); x})
                e2 <- elsa(o1,d1=d1,d2=d2,nc=nc,categorical=categorical,dif=dif,cells=cells)
                ee <- e1 - e2
                ee <- ifelse(ee > 0,1,0)
                o2 <- o2 + ee
              }
              rm(e1,e2,ee,o1)
              o2 <- (o2+1) / (n+1)
            }
            o2
          }
)

