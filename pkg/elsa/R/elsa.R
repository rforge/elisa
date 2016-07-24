# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2016
# Version 1.3
# Licence GPL v3 


.Filter<-function(r,d1=0,d2) {
  c<- d2%/%r
  x<- y<- seq(-c,c,1)
  eg<- expand.grid(x,y)
  eg[1]<- -eg[1]
  eg[3]<- sqrt(eg[1]^2+eg[2]^2)
  ndim<- c*2+1
  m<-matrix(eg[,3],ncol=ndim,nrow=ndim)*r
  mw<- which(m > d1 & m <= d2)
  m[mw] = 1 
  m[-mw] = NA
  mw <- trunc(length(m)/2)+1
  m[mw]<- 1
  
  w <- ncol(m)
  mr <- matrix(ncol=w,nrow=w) 
  tr <- trunc(w/2)
  
  for (i in 1:w) {
    mr[i,] <-tr 
    tr <- tr - 1
  }
  
  rc = cbind(r=as.vector(m * mr),c=as.vector(m * t(mr)))
  rc=rc[apply(rc,1,function(x) (all(!is.na(x)))),]
  
  return (list(w,rc))
}

#----------
.is.categorical <- function(x) {
  o <- unique(x)
  if (length(o) == length(round(o)) && length(o) < 100) o <- TRUE
  else o <- FALSE
  o
}
#----------
.checkDif <- function(dif,classes) {
  nc <- length(classes)
  classes <- as.character(classes)
  if(is.list(dif)) {
    if(!is.null(names(dif))) {
      if (!all(unlist(lapply(classes,function(x) {x %in% names(dif)})))) stop("categories in the dif argument do not match with the categories in the input layer!")
      dif <- dif[classes]
    } else {
      if (length(dif) == nc) names(dif) <- classes
      else stop("dif content does not match with the categories in the input layer!")
    }
    if (!all(unlist(lapply(classes,function(x) {length(dif[[x]] == nc)})))) stop("the length of categories' differences in the dif argument does not match with the number of categories!")
    dif <- as.vector(as.matrix(as.data.frame(dif)))
  } else if (is.matrix(dif) || is.data.frame(dif)) {
    if (ncol(dif) != nrow(dif) || ncol != nc) stop("number of rows and columns in dif should be the same, equal to the number of categories!")
    if (!all(unlist(lapply(colnames(dif),function(x) {x %in% row.names(dif)})))) {
      if (!all(unlist(lapply(colnames(dif),function(x) {x %in% classes})))) {
        dif <- as.vector(as.matrix(as.data.frame(dif)))
      } else {
        row.names(dif) <- colnames(dif)
        dif <- dif[classes,classes]
        dif <- as.vector(as.matrix(as.data.frame(dif)))
      }
    } else {
      if (!all(unlist(lapply(colnames(dif),function(x) {x %in% classes})))) stop("colnames of dif do not match with the name of categories!")
      else {
        dif <- dif[classes,classes]
        dif <- as.vector(as.matrix(as.data.frame(dif)))
      }
    }
    
  } else if (is.vector(dif)) {
    if (length(dif) != nc*nc) stop("dif argument does not have an appropriate structure!")
  } else stop("dif argument does not have an appropriate structure!")
  return(dif)
}
#----------
if (!isGeneric("elsa")) {
  setGeneric("elsa", function(x,d,nc,categorical,dif,...)
    standardGeneric("elsa"))
}


setMethod('elsa', signature(x='RasterLayer'), 
          function(x,d,nc,categorical,dif,cells,filename,...) {
            if (missing(d)) d <- res(x)[1]
            if (missing(filename)) filename=''
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
            w <-.Filter(r=res(x)[1],d1=0,d2=d)
            fdim <- w[[1]]
            w <- w[[2]]
            
            if (fdim < 3) stop("d must be at least equal to the input raster resolution!")
            
            if (!categorical) x <- categorize(x,nc)
            
            out <- raster(x)
            ncl <- ncol(out)
            nrw <- nrow(out)
            filename=trim(filename)
            
            if (canProcessInMemory(out)) {
              if (categorical) {
                if (missing(cells)) {
                  out[] <- .Call('elsac', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif)
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  out <- .Call('elsac_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cells))
                }
              } else {
                if (missing(cells)) {
                  out[] <- .Call('elsa', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]))
                  if (filename != '') out <- writeRaster(out, filename, ...)
                } else {
                  out <- .Call('elsa_cell', as.integer(x[]), as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(cells))
                }
              }
            } else {
              tr <- blockSize(out, minblocks=3, minrows=fdim)
              pb <- pbCreate(tr$n, label='ELSA',...)
              addr <- floor(fdim / 2)
              
              if (missing(cells)) {
                out <- writeStart(out, filename)
                v <- as.integer(getValues(x, row=1, nrows=tr$nrows[1]+addr))
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]))
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif)
                }
                ex <- length(v) - (addr * ncl)
                out <- writeValues(out, v[1:ex], 1)
                
                for (i in 2:(tr$n-1)) {
                  v <- as.integer(getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr)))
                  if (!categorical) {
                    v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]))
                  } else {
                    v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif)
                  }
                  st <- (addr * ncl)+1
                  ex <- length(v) - (addr * ncl)
                  out <- writeValues(out, v[st:ex], tr$row[i])
                  pbStep(pb)
                }
                
                i <- tr$n
                v <- as.integer(getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr))
                if (!categorical) {
                  v <- .Call('elsa', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]))
                } else {
                  v <- .Call('elsac', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif)
                }
                st <- (addr * nc)+1
                ex <- length(v)
                out <- writeValues(out, v[st:ex], tr$row[i])
                pbStep(pb)
                out <- writeStop(out)      
                pbClose(pb)  
              } else {
                v <- as.integer(getValues(x, row=1, nrows=tr$nrows[1]+addr))
                cls <- cells[which(cells <= (tr$nrows[1]) * ncl)]
                if (length(cls) > 0) {
                  if (!categorical) {
                    v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls))
                  } else {
                    v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls))
                  }
                  out <- c(out, v)
                }
                
                for (i in 2:(tr$n-1)) {
                  cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                  if (length(cls) > 0) {
                    cls <- cls - ((tr$row[i]-addr-1)*ncl)
                    v <- as.integer(getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+(2*addr)))
                    if (!categorical) {
                      v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls))
                    } else {
                      v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls))
                    }
                    out <- c(out, v)
                    pbStep(pb)
                  }
                }
                
                i <- tr$n
                cls <- cells[which(cells > ((tr$nrow[i] - 1) * ncl) & cells <= ((tr$nrows[i]+ tr$nrows[i] - 1) * ncl))]
                cls <- cls - ((tr$row[i]-addr-1)*ncl)
                if (length(cls) > 0) {
                  v <- as.integer(getValues(x, row=tr$row[i]-addr, nrows=tr$nrows[i]+addr))
                  if (!categorical) {
                    v <- .Call('elsa_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]), as.integer(cls))
                  } else {
                    v <- .Call('elsac_cell', v, as.integer(ncl), as.integer(nrw), as.integer(nc), as.integer(w[,1]), as.integer(w[,2]),as.integer(classes),dif, as.integer(cls))
                  }
                  out <- c(out, v)
                  pbStep(pb)
                  pbClose(pb)  
                }
              }
            }
            return(out)
          }
)  

#---------------
# 
# setMethod('elsa', signature(x='SpatialPointsDataFrame'), 
#           function(x,d,nc,categorical,dif,zcol,...) {
#             if (missing(d)) stop('d is missed!')
#             if (missing(zcol)) {
#               if (ncol(x@data) > 1) stop("zcol should be specified!")
#               else zcol <- 1
#             } else if (is.character(zcol)) {
#               w <- which(colnames(x@data) == zcol[1])
#               if (w != 1) stop('the specified variable in zcol does not exist in the data')
#               zcol <- w
#             } else if (is.numeric(zcol)) {
#               zcol <- zcol[1]
#               if (zcol > ncol(x@data)) stop('the zcol number is greater than the number of columns in data!')
#             } else stop("zcol should be a character or a number!")
#             
#             xy <- coordinates(x)
#             yy <- xy[,2]
#             xx <- xy[,1]
#             x <- x@data[,zcol]
#             rm(xy)
#             
#             if (is.character(x) || is.factor(x)) {
#               x <- as.character(x)
#               if (!missing(categorical) && !categorical) warning("you specified a categorical variable, so categorical changed to TRUE!")
#               categorical <- TRUE
#             }
#             
#             if (!missing(nc)) {
#               if (missing(categorical)) {
#                 if (missing(dif)) categorical <- FALSE
#                 else {
#                   categorical <- TRUE
#                   cat("input data is considered categorical, and nc is ignored!\n")
#                 }
#               } 
#             } else {
#               if (missing(categorical) && !missing(dif)) categorical <- TRUE
#             }
#             #----
#             if (missing(categorical) || !is.logical(categorical)) {
#               # guessing whether the layer is categorical:
#               if (.is.categorical(x)) {
#                 categorical <- TRUE
#                 cat("the specified variable is considered as categorical...\n")
#               } else {
#                 categorical <- FALSE
#                 cat("the specified variable is considered continuous...\n")
#               }
#             }
#             #----
#             if (!categorical && missing(nc)) {
#               nc <- nclass(x)
#             } else if (categorical) {
#               classes <- unique(x)
#               nc <- length(classes)
#             }
#             #-----
#             
#             if (categorical) {
#               if (missing(dif)) {
#                 dif <- rep(1,nc*nc)
#                 for (i in 1:nc) dif[(i-1)*nc+i] <-0
#               } else {
#                 dif <- .checkDif(dif,classes)
#               }
#             }
#             #-----
#             
#             if (!categorical) x <- categorize(x,nc)
#             
#             if (categorical) {
#               
#             } else {
#               .Call('elsa_vector', as.integer(x), as.vector(xx), as.vector(yy), as.integer(nc), as.numeric(d))
#             }
#             
#           }
# )  




setMethod('elsa', signature(x='SpatialPointsDataFrame'), 
          function(x,d,nc,categorical,dif,zcol,...) {
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, d[1])
            d <- d@neighbours
            
            if (missing(zcol)) {
              if (ncol(x@data) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(colnames(x@data) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x@data)) stop('the zcol number is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            
            x <- x@data[,zcol]
            
            if (is.character(x) || is.factor(x)) {
              x <- as.character(x)
              if (!missing(categorical) && !categorical) warning("you specified a categorical variable, so categorical changed to TRUE!")
              categorical <- TRUE
            }
            
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
                cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                cat("the specified variable is considered continuous...\n")
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
            
            if (!categorical) x <- categorize(x,nc)
            
            
            if (categorical) {
              .Call('elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif)
            } else {
              .Call('elsa_vector', as.integer(x), d, as.integer(nc))
            }
            
          }
)  


setMethod('elsa', signature(x='SpatialPolygonsDataFrame'), 
          function(x,d,nc,categorical,dif,zcol,method,...) {
            if (missing(d)) stop('d is missed!')
            else if (!class(d) %in% c('numeric','integer','neighbours')) stop('d should be either a number (distance) or an object of class neighbours (created by dneigh function')
            
            if (missing(method)) method <- 'centroid'
            
            if (!inherits(d,'neighbours')) d <- dneigh(x, d[1],method = method)
            d <- d@neighbours
            
            
              if (missing(zcol)) {
              if (ncol(x@data) > 1) stop("zcol should be specified!")
              else zcol <- 1
            } else if (is.character(zcol)) {
              w <- which(colnames(x@data) == zcol[1])
              if (w == 0) stop('the specified variable in zcol does not exist in the data')
              zcol <- w
            } else if (is.numeric(zcol)) {
              zcol <- zcol[1]
              if (zcol > ncol(x@data)) stop('the zcol number is greater than the number of columns in data!')
            } else stop("zcol should be a character or a number!")
            
            x <- x@data[,zcol]
            
            if (is.character(x) || is.factor(x)) {
              x <- as.character(x)
              if (!missing(categorical) && !categorical) warning("you specified a categorical variable, so categorical changed to TRUE!")
              categorical <- TRUE
            }
            
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
                cat("the specified variable is considered as categorical...\n")
              } else {
                categorical <- FALSE
                cat("the specified variable is considered continuous...\n")
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
            
            if (!categorical) x <- categorize(x,nc)
            
            if (categorical) {
              .Call('elsac_vector', as.integer(x), d, as.integer(nc), as.integer(classes),dif)
            } else {
              .Call('elsa_vector', as.integer(x), d, as.integer(nc))
            }
            
          }
)  
