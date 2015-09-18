# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

#' categorizing continious data
#' 
#' A function to categorize numerical data in the form of a raster layer, or a vector
#' 
#' @aliases categorize, categorize,Raster-method, categorize,numeric-method
#' 
#' @usage categorize(x,nc,filename,...)
#' 
#' @param x A \code{RasterLayer} or a \code{numerical vector}
#' @param nc Number of classes, if missing, it will be automatically detected
#' @param filename A character specifying the output Raster file name (if x is RasterLayer)
#' @param ... Additional parameters for \code{\link[raster]{writeRaster}} function (if filename is specified)
#' 
#' @return An object with the same class as \code{x}
#' @author Babak Naimi \email{naimi.b@@gmail.com}
#' 
#' @examples 
#' categorize(1:10,3)
#' file <- system.file('external/dem_example.grd',package='elsa')
#' r <- raster(file)
#' plot(r,main='a continuous raster layer')
#' rc <- categorize(r,nc=4)
#' plot(rc, main=') 



                        
if (!isGeneric("categorize")) {
  setGeneric("categorize", function(x, ...)
    standardGeneric("categorize"))
}	

setMethod('categorize', signature(x='RasterLayer'), 
          function(x,nc,filename='',...)  {
            if (missing(nc)) {
              nc <- nclass(x)
              cat(paste("the optimum number of class has been identified as ",nc,"!\n"))
            }
            
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- cellStats(x,'range')
              n <- (r[2] - r[1])/ nc
              nc <- seq(r[1],r[2],n)
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            out <- raster(x)
            #-----
            if (canProcessInMemory(out)) {
              out[] <- .Call('categorize', as.vector(x[]), as.vector(nc))
              if (filename != '') out <- writeRaster(out, filename, ...)
            } else {
              out <- writeStart(out, filename,...)
              tr <- blockSize(out, minblocks=3, minrows=fdim)
              pb <- pbCreate(tr$n, label='categorize',...)
              
              for (i in 1:tr$n) {
                v <- getValues(x, row=tr$row[i], nrows=tr$nrows[i])
                v <- .Call('categorize', v, as.vector(nc))
                out <- writeValues(out, v, 1)
                pbStep(pb)
              }
              out <- writeStop(out)      
              pbClose(pb)
            }
            return(out)
          }
)
            



setMethod('categorize', signature(x='numeric'), 
          function(x,nc)  {
            if (missing(nc)) {
              stop("number of classes or a verctor including the break values should be specified...!")
            }
            
            if (length(nc) == 1) {
              if (nc < 2) stop("nclass should be 2 or greater!")
              r <- range(x,na.rm=TRUE)
              n <- (r[2] - r[1])/ nc
              nc <- seq(r[1],r[2],n)
              nc[1] <- nc[1] - n
              if (nc[length(nc)] < r[2]) nc[length(nc)] <- r[2]
            }
            .Call('categorize', as.vector(x), as.vector(nc))
          }
)

