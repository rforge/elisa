\name{entrogram}

\alias{entrogram}
\alias{entrogram,RasterLayer-method}


\title{Empirical Entrogram for Spatial Data}

\description{
Compute sample (empirical) variogram from spatial data. The function returns a binned entrogram and an entrogram cloud.
}

\usage{
entrogram(x, width, cutoff,...) 
}

\arguments{
  \item{x}{ a spatial object (\code{RasterLayer} or \code{SpatialPoints} or \code{SpatialPolygons}) }
  \item{width}{ the lag size (width of subsequent distance intervals) into which cell pairs are grouped for semivariance estimates. If missing, the cell size (raster resolution) is assigned. }
  \item{cutoff}{ spatial separation distance up to which cell pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three. }
  \item{...}{ }
  
}
\details{
Variograms are widely used for exploring spatial structure in a single variable. Formally, it is defined as half the expected squared difference (half the variance of the difference) in the variable value at a specific geographical separation. A variogram summarizes the spatial relations in the data, and can be used to understand within what range (distance) the data is spatially autocorrelated.
}

\value{
\item{Entrogram}{}
}

\references{
Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

}

\keyword{spatial}
