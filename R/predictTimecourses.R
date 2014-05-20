predictTimecourses <- function(res.fits){
	
	fits = res.fits$fitted
	timepoints = res.fits$time
	nc = res.fits$cols.mat
	n.phi = length(res.fits$phi)
	
	lower.standard = timepoints[1]
	if(lower.standard == 0) lower.standard = 1
	upper.standard= timepoints[length(timepoints)-1]
	
	phi.resolution = res.fits$phi[2]-res.fits$phi[1]
	time.resolution = timepoints[length(timepoints)]/length(timepoints)
	if(phi.resolution < time.resolution){timepoints=seq(timepoints[1],timepoints[length(timepoints)],by=phi.resolution)}
	
	res = list()
	ids = c()
	for(i in 1:length(fits)){
		f = fits[[i]]
		ids = c(ids,f$ID)
		if(is.null(f$psi)){
			wave = makewave(f,timepoints,lower=lower.standard,upper=upper.standard,fun=fun.lognormal)
		} else {
			wave = makewave(f,timepoints,lower=lower.standard,upper=upper.standard,psi=f$psi,fun=fun.lognormal.psi)
		}
		wave = wave*f$a+f$b
		res[[length(res)+1]] = wave	
	}

	names(res) = ids
	res = t(as.data.frame(res))
	colnames(res) = timepoints
	if(dim(res)[1] == 1) {return(as.numeric(res))}
	else { return(res)}
}
