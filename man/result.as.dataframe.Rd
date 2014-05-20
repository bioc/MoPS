\name{result.as.dataframe}
\alias{result.as.dataframe}
\title{
Parameters of fitted periodic time courses.
}
\description{
Converts the result list derived by MoPS function fit.periodic() to a data.frame. 
}
\usage{
result.as.dataframe(result.list)
}
\arguments{
  \item{result.list}{
List of best fitting parameters returned by fit.periodic().
}
}
\details{
This function takes as input the result list from MoPS function fit.periodic() and extracts the time course specific optimal parameters.

}
\value{
data.frame containing the best fitting periodic parameters for each time series (rows):

columns:
ID : unique identifier
score : log-likelihood for periodic behaviour
phi : phase
lambda : period length
sigma : attenuation of the signal along the complete time series
mean : mean 
amplitude : amplitude

}
\author{
Philipp Eser, Achim Tresch
}

\examples{

y = 2*sin(seq(0,6*pi,length.out=50)+rnorm(50))
res = fit.periodic(y)
result.as.dataframe(res)

}

