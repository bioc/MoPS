\name{predictTimecourses}
\alias{predictTimecourses}
\title{
Prediction of periodic time courses. 
}
\description{
Function that predicts periodic time courses using parameters identified by fit.periodic(). 
}
\usage{
predictTimecourses(res.fits)
}
\arguments{
  \item{res.fits}{
List object returned by fit.periodic().
}
  
}
\details{
This function takes as input the result list from MoPS function fit.periodic() and creates a list of best fitting time courses.
The input list also contains information about the screening parameters, which is used in the generation of predicted time courses.
}
\value{
a numeric matrix containing the predicted values. 
The number of rows equals the number of rows of the original data matrix, the number of columns equals the number of screened phases. 
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

