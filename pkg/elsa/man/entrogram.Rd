\name{entrogram}

\alias{entrogram}
\alias{entrogram,RasterLayer-method}


\title{Empirical Entrogram for Raster Data}

\description{
Compute sample (empirical) variogram from raster data. The function returns a binned variogram and a variogram cloud.
}

\usage{
entrogram(x, width, cutoff,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer}) }
  \item{width}{ the lag size (width of subsequent distance intervals) into which cell pairs are grouped for semivariance estimates. If missing, the cell size (raster resolution) is assigned. }
  \item{cutoff}{ spatial separation distance up to which cell pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three. }
  \item{...}{ }
  
}
\details{
Variograms are widely used for exploring spatial structure in a single variable. Formally, it is defined as half the expected squared difference (half the variance of the difference) in the variable value at a specific geographical separation. A variogram summarizes the spatial relations in the data, and can be used to understand within what range (distance) the data is spatially autocorrelated. Naimi et al. (2011) linked this range to the impact of positional uncertainty on the performance of species distribution models (SDMs). Based on that study, examining variogram to find the effective autocorrelation range in predictors gives insight into whether predictions by SDMs are likely to be affected by the uncertainty in the sample locations (see Naimi et al. 2011, for more information).
}

\value{
\item{RasterVariogram}{}
}

\references{
Naimi, B., Skidmore, A.K, Groen, T.A., Hamm, N.A.S. 2011. Spatial autocorrelation in predictors reduces the impact of positional uncertainty in occurrence data on species
distribution modelling, Journal of biogeography. 38: 1497-1509.
}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

}

\keyword{spatial}
