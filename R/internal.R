


createLinearTimeCourses <- function(timepoints){
	
	npoints = length(timepoints)
	
	uplinear = seq(from=-1,to=1,length=npoints)
	uplinear = uplinear-mean(uplinear)
	uplinear = uplinear/sqrt(sum(uplinear^2))
	
	m = diag(npoints)
		
	fm = cbind(rep(1/sqrt(npoints),npoints),uplinear,rev(uplinear),m,-m)
	return(fm)
}

decimpart <- function(n){n-floor(n/(2*pi))*2*pi}


extend <- function(v,extendgrid){ 
	extendable = extendgrid[max(v) <= extendgrid] 
	mat = rbind(matrix(rep(v,length(extendable)),ncol=length(extendable)),extendable)
	return(mat)
}

fun.lognormal.psi <- function(x,lambda,phi,sigma,t=0,psis=psi){
	
	logsigma = log(exp(log(sigma^2)-2*log(lambda))+1)
	logmu = log(lambda)-logsigma/2
	psi = pathfunction(psis)
	#convert phase to relative circle measure
	phi = 2*pi*phi/lambda
	res = cos(psi(decimpart(2*pi*t/x-phi)))*dlnorm(x,meanlog=logmu,sdlog=sqrt(logsigma))
	
	return(res)
}

fun.lognormal <- function(x,lambda,phi,sigma,t=0){
	logsigma = log(exp(log(sigma^2)-2*log(lambda))+1)
	logmu = log(lambda)-logsigma/2
	return( cos((2*pi*(t-phi))/x)*dlnorm(x,meanlog=logmu,sdlog=sqrt(logsigma)) )
}

makewave <-	function(x,times,lambda=x[["lambda"]],phi=x[["phi"]],sigma=x[["sigma"]],
				lower,upper,fun=fun.lognormal,psi=NULL){
	wave = numeric(length(times))
	
	if(is.null(psi)){
		if(sigma == 0){
			wave = cos((2*pi*(times-phi))/lambda)
		} else {
			for (k in 1:length(times)){
				wave[k] = integrate(fun, lower = lower,upper = upper,
						lambda = lambda, phi = phi, sigma = sigma,t = times[k])$value
			}
		}
	} else {
		if(sigma == 0){
			psi = pathfunction(psi)
			wave = cos(psi(decimpart(2*pi*times/lambda-phi)))
			
		} else {
			for (k in 1:length(times)){
				wave[k] = integrate(fun, lower = lower,upper = upper,
						lambda = lambda, phi = phi, sigma = sigma,psis=psi,t = times[k],stop.on.error=FALSE)$value
			}
		}
	}
	return(wave)
}


pathfunction <-	function(y,x=NULL){
	if (is.null(x)) x = seq(0,2*pi,length=length(y)+2)[-c(1,length(y)+2)]
	return( approxfun(c(0,x,2*pi),c(0,y,2*pi), method = "linear") )
}

