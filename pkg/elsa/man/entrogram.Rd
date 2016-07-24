\name{entrogram}

\alias{entrogram}
\alias{entrogram,RasterLayer-method}


\title{Entrogram for Spatial Data}

\description{
Compute sample (empirical) entrogram from spatial data. The function returns a binned entrogram and an entrogram cloud.
}

\usage{
entrogram(x, width, cutoff,...) 
}

\arguments{
  \item{x}{ a spatial object (\code{RasterLayer} or \code{SpatialPoints} or \code{SpatialPolygons}) }
  \item{width}{ the lag size (width of subsequent distance intervals) into which cell pairs are grouped for semivariance estimates. If missing, the cell size (raster resolution) is assigned. }
  \item{cutoff}{ spatial separation distance up to which cell pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three. }
  \item{...}{Additional arguments passed to \code{\link[elsa]{elsa}} function }
  
}
\details{
Entrogram is a variogram-like graph to explore spatial structure in a single variable. An entrogram summarizes the spatial relations in the data, and can be used to understand within what range (distance) the data is spatially autocorrelated.
}

\value{
\item{Entrogram}{an object containing entrogram cloud and the entrogram within each distance interval}
}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{
\dontrun{
file <- system.file('external/dem_example.grd',package='elsa')
r <- raster(file)
plot(r,main='a continuous raster map')

en <- entrogram(r, width=2000)

plot(en)
}


}

\keyword{spatial}
