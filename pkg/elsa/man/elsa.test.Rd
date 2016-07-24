\name{elsa.test}

\alias{elsa.test}
\alias{elsa.test,RasterLayer-method}



\title{Elsa test for local spatial autocorrelation}

\description{
This function can be used to convert a data.frame including the values specifying the differences (digree of dissimilarities) between the classes (categories) in a categoricl map, to a list.
}

\usage{
elsa.test(x, d, n, nc, categorical, dif,...)
}

\arguments{
  \item{x}{ a data.frame containing digrees of dissimilarity between different categories}
  \item{d}{the local distance}
  \item{n}{number of simulation, default is 99}
  \item{nc}{number of classes (only if x is a continuous variable); if not specified, it is estimated using nclass function}
  \item{categorical}{logical, specifies whether x is a categorical; if not specified, it is guessed by the function}
  \item{dif}{the level of dissimilarities between different categories (only if x is a categorical variable); see \code{\link[elsa]{dif2list}} for more details}
  \item{...}{Aditional arguments passed to elsa function}
}
\details{
This function test how significant the local spatial autocorrelation is at each location, so it generates a p-value at each location through a Monte Carlo simulation and a non-parametric approach. See the reference for the details about the method.

}

\value{
An object same as the input (\code{x})

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

et <- elsa.test(r,d=2000,n=99)

plot(et)

}

}

\keyword{spatial}
