# Author: Babak Naimi, naimi.b@gmail.com
# Date :  August 2014
# Version 1.0
# Licence GPL v3 

# based on the functions poly2nb and dnearneigh in spdep package (Roger Bivand):

.qintersect <- function (x, y) {
  as.integer(y[match(x, y, 0L)])
}


.findInbox <- function (i, sp, bigger = TRUE) {
  n <- dim(sp$bb)[1]
  tmp <- vector(mode = "list", length = 4)
  tmp[[1]] <- sp$rbxv[sp$mbxv[i]:(n * 2)]
  tmp[[1]] <- tmp[[1]][which(tmp[[1]] > n)] - n
  tmp[[2]] <- sp$rbyv[sp$mbyv[i]:(n * 2)]
  tmp[[2]] <- tmp[[2]][which(tmp[[2]] > n)] - n
  tmp[[3]] <- sp$rbxv[1:sp$mbxv[i + n]]
  tmp[[3]] <- tmp[[3]][which(tmp[[3]] <= n)]
  tmp[[4]] <- sp$rbyv[1:sp$mbyv[i + n]]
  tmp[[4]] <- tmp[[4]][which(tmp[[4]] <= n)]
  lentmp <- order(sapply(tmp, length))
  result <- .qintersect(tmp[[lentmp[2]]], tmp[[lentmp[1]]])
  result <- .qintersect(tmp[[lentmp[3]]], result)
  result <- .qintersect(tmp[[lentmp[4]]], result)
  if (bigger) {
    result <- result[which(result > i)]
  }
  return(sort(result))
}

.dnn.poly <- function(x,queen=TRUE) {
  n <- length(slot(x, "polygons"))
  if (n < 1) stop("non-positive number of entities")
  regid <- row.names(x)
  if (is.null(regid)) regid <- as.character(1:n)
  
  xpl <- slot(x, "polygons")
  xxpl <- vector(mode = "list", length = length(xpl))
  
  for (i in 1:length(xpl)) {
    xpli <- slot(xpl[[i]], "Polygons")
    zz <- lapply(xpli, function(j) slot(j, "coords")[-1, ])
    xxpl[[i]] <- do.call("rbind", zz)
  }
  nrs <- sapply(xxpl, nrow)
  vbsnap <- c(-d, d)
  dsnap <- as.double(d)
  bb <- t(sapply(xxpl, function(x) {
    rx <- range(x[, 1]) + vbsnap
    ry <- range(x[, 2]) + vbsnap
    c(rbind(rx, ry))
  }))
  genBBIndex <- function(bb) {
    n <- nrow(bb)
    bxv <- as.vector(bb[, c(1, 3)])
    byv <- as.vector(bb[, c(2, 4)])
    obxv <- order(bxv)
    rbxv <- c(1:(n * 2))[obxv]
    mbxv <- match(1:(n * 2), obxv)
    obyv <- order(byv)
    rbyv <- c(1:(n * 2))[obyv]
    mbyv <- match(1:(n * 2), obyv)
    return(list(bb = bb, bxv = bxv, byv = byv, obxv = obxv,obyv = obyv, mbxv = mbxv, mbyv = mbyv, rbyv = rbyv, rbxv = rbxv))
  }
  BBindex <- genBBIndex(bb)
  foundInBox <- lapply(1:(n - 1), function(i) .findInbox(i, BBindex))
  nfIBB <- sum(sapply(foundInBox, length))
  criterion <- ifelse(queen, 0, 1)
  ans <- .Call("poly_loop2", as.integer(n), foundInBox, bb, xxpl, as.integer(nrs), as.double(dsnap), as.integer(criterion), as.integer(nfIBB), PACKAGE = "elsa")
  ans
}
#------
.dnn.xy <- function(x,d,longlat) {
  .Call("dnn", x[,1], x[,2], d, as.integer(longlat), PACKAGE = "elsa")
}

######################################
.is.projected <- function(x) {
  if (inherits(x,'Spatial')) {
    if (!is.na(is.projected(x))) {
      is.projected(x)
    } else {
      all(bbox(x)[1,] <= 180) & all(bbox(x)[1,] >= -180) & all(bbox(x)[2,] <= 90) & all(bbox(x)[2,] >= -90)
    }
  } else if (inherits(x,'matrix') || inherits(x,'data.frame')) {
    all(range(x[,1],na.rm=TRUE) <= 180) & all(range(x[,1],na.rm=TRUE) >= -180) & all(range(x[,2],na.rm=TRUE) <= 90) & all(range(x[,2],na.rm=TRUE) >= -90)
  } 
}


if (!isGeneric("dneigh")) {
  setGeneric("dneigh", function(x, d,longlat,method,...)
    standardGeneric("dneigh"))
}



setMethod('dneigh', signature(x='SpatialPoints',d='numeric'), 
          function(x, d, longlat,...) {
            if (missing(longlat) || is.null(longlat) || !is.logical(longlat)) longlat <- .is.projected(x)
            x <- coordinates(x)
            if (nrow(x) < 1) stop("no records in x")
            if (ncol(x) > 2) stop("Only 2D data accepted")
            z <- .dnn.xy(x,d,longlat)
            if (all(sapply(z,is.null))) stop('There is no links within the specified distance!')
            z <- new('neighbours',distance1=0,distance2=d,neighbours=z)
            return(z)
          }
)



setMethod('dneigh', signature(x='SpatialPolygons',d='numeric'), 
          function(x, d, longlat,method,...) {
            if (missing(longlat) || is.null(longlat) || !is.logical(longlat)) longlat <- .is.projected(x)
            
            if (missing(method) || is.null(method)) method <- 'bound'
            else {
              if (tolower(method)[0] %in% c('bnd','bound','boundary','bond','b')) method <- 'bound'
              else if (tolower(method)[0] %in% c('center','centre','cent','cnt','c','centroid','centriod','ce','cen')) method <- 'centroid'
              else {
                warning('method is not recognized; the default (bound) is used!')
                method <- 'bound'
              }
            }
            
            if (method == 'bound') {
              dot <- list(...)
              if ('queen' %in% names(dot) && is.logical(dot[['queen']])) queen <- dot[['queen']]
              else queen <- TRUE
              z <- .dnn.poly(x,d=d,queen=queen)
              attributes(z) <- NULL
            } else {
              x <- coordinates(x)
              if (nrow(x) < 1) stop("no records in x")
              z <- .dnn.xy(x,d,longlat)
            }
            if (all(sapply(z,is.null))) stop('There is no links within the specified distance!')
            z <- new('neighbours',distance1=0,distance2=d,neighbours=z)
            return(z)
          }
)

setMethod('dneigh', signature(x='data.frameORmatrix',d='numeric'), 
          function(x, d, longlat,...) {
            if (nrow(x) < 1) stop("no records in x")
            if (missing(longlat) || is.null(longlat) || !is.logical(longlat)) longlat <- .is.projected(x)
            z <- .dnn.xy(as.matrix(x),d,longlat)
            if (all(sapply(z,is.null))) stop('There is no links within the specified distance!')
            z <- new('neighbours',distance1=0,distance2=d,neighbours=z)
            return(z)
          }
)
