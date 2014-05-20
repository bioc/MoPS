\name{fit.periodic}
\alias{fit.periodic}
\title{
Fitting of periodic curves to time series data
}
\description{
Function that fits periodic curves to each time course in a numeric matrix or ExpressionSet object. 
It determines the best fitting time courses from an exhaustive set of periodic and non-periodic test functions.
These are created either automatically or with a user-specified parameter set (see Arguments).
Fitting is done with standard linear regression. 

In addition to the best fitting periodic curve, the best-fitting linear time course is estimated.
A periodicity score is derived based on the godness-of-fit ratio between periodic and non-periodic fits.

The function returns a list containing the best fitting parameters for each time series. The function result.as.dataframe() converts the result into a data.frame.  
}
\usage{
fit.periodic(mat,timepoints=NULL,phi=NULL,lambda=NULL,sigma=NULL,psi=NULL,weights=NULL)
}
\arguments{
  \item{mat}{
a numeric matrix containing individual measurements (rows) across a time series (columns). mat can also supplied as an ExpressionSet object.  
}
  \item{timepoints}{
optional numeric vector corresponding to measurement timepoints. If NULL, timepoints are initialized as 1:ncol(mat).
}
  \item{phi}{
optional numeric vector specifying all possible phases of periodic test functions (phase = time where periodic curve is maximal).  
}
  \item{lambda}{
optional numeric vector specifying all possible period lengths of periodic test functions.
}
  \item{sigma}{
optional numeric vector specifying the magnitude of dampening of the signal along the time course.
}
 \item{psi}{
optional positive integer defining the level of flexibility to shape the test functions. Recommended values for psi are 3 or 4. psi > 4 results in a tremendous increase in runtime.
}
 \item{weights}{
optional numeric matrix of weights to be used in the fitting. If non-NULL, weighted least squares is used, otherwise ordinary least squares is used.
}
}
\details{
The input data needs to be a numeric matrix containing in each row a time series of measurements.

The function can take an optional numeric matrix of weights as input that is used in the fitting process.
This matrix needs to have the same dimensions as the input data matrix.
If weights are supplied, weighted least squares is used otherwise ordinary least squares is used.
This option is useful if the size of the measurement error is not constant for all measurements.

Note that this function uses all possible parameter combinations to create periodic test functions.
This can be very time consuming if the user chooses wide parameter ranges as input.
If possible, the user should specify meaningful ranges with a moderate spacing between values (see also the MoPS vignette).      
}
\value{

fit.periodic() returns a list object containing information about
the fitting results for each input time series and the parameter ranges used in the screening.

The first slot of the result object contains the following values for each time series: 
\item{$ID}{unique id}
\item{$is.wPeriodic}{TRUE if $minLossPeriodic < $minLossNonPeriodic }
\item{$minLossPeriodic}{loss of best periodic fit}
\item{$minLossNonPeriodic}{loss of best non-periodic fit}
\item{$phi}{phase} 
\item{$psi}{variable sampling points of best fitting psi transformation}
\item{$lambda}{period length}
\item{$sigma}{signal attenuation along the time series}
\item{$a.coef}{coefficient a from linear model (amplitude)}
\item{$b.coef}{coefficient b from linear model (mean)}
 
The remaining slots contain the following values: 
\item{$time}{measurement time points}
\item{$cols.mat}{number of columns of the input data matrix}
\item{$phi}{all screened phi values}
\item{$lambda}{all screened lambda values}
\item{$sigma}{all screened sigma values}
 
For convenient sorting or filtering, this list can be converted to a data.frame with the function result.as.dataframe(). 
}
\references{
Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.  
}
\author{
Philipp Eser, Achim Tresch
}
\examples{

x = seq(0,40,by=1) # time points

## create 10 periodic time series with added noise
mat.p = matrix(rep(x,10),nrow=10,ncol=length(x),byrow=TRUE)
y = -seq(1:10)
mat.p = apply(mat.p,2,function(x){
	y = sin(pi*(x/41*6)+y)+rnorm(length(x),sd=1)
	})

## add 10 non-periodic noisy time series
mat.nonP = matrix(rep(x,10),nrow=10,ncol=length(x),byrow=TRUE)
mat.nonP = apply(mat.nonP,2,function(x){
	y = rnorm(length(x),sd=1)
	})
	
mat = rbind(mat.p,mat.nonP)

res = fit.periodic(mat,phi=seq(0,20,1),lambda=seq(1,20,1))
time.courses = predictTimecourses(res)

plot(mat[1,],type="l",main="",xlab="",ylab="")
points(time.courses[1,],type="l",col="limegreen",lwd=2)

}