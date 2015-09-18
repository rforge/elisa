\name{elsa}

\alias{elsa}
\alias{elsa,Raster,missing-method}
\alias{elsa,Raster,SpatialPoints-method}
\alias{lisa,Raster,SpatialPointsDataFrame-method}

\title{Entropy-based Local indicators of Spatial Sssociation}

\description{
Calculate ELSA indicator for either a categorical or continuous spatial dataset.
}

\usage{
elsa(x,d,nc,categorical,dif,...) 
}

\arguments{
  \item{x}{ a raster object (\code{RasterLayer}}
  \item{d}{ local distance  }
  \item{nc}{ optional, for continuous data it specifies the number of classes through categorizing the variable. If missed, it is automatically calculated (recommended)  }
  \item{categorical}{ logical, specified whether x is a continuous or categorical. If missed the function tries to detect it}
  \item{dif}{the difference between categories, only for categorical }
  \item{...}{ additional arguments including \code{cells} (when the elsa is needed for only some cells in the raster data), and \code{filename} (the name of the raster file when the output should be written) }
}
\details{
{dif} can be used when categorical values are sorted into hierarchical system (e.g., CORINE land cover). This make it possible to difine different weights of similarity between each pairs of categories when the level of similarity is not the same between different classes in the variable. For example, two categories belong to two forest types are more similar than two categories, one a forest type and the other one an agriculture type. So, it can take this differences into account when the spatial autocorrelation for categorical variables is quantified.
}

\value{
\item{RasterLayer}{if \code{x} is a \code{RasterLayer}

}

\references{Naimi. B. et al. (under review) ELSA: An Entropy-based Local indicators of Spatial Association, Geographical Analysis;

}

\author{Babak Naimi \email{naimi.b@gmail.com}

\url{http://r-gis.net}
}


\examples{

}

\keyword{spatial}
