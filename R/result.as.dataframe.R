result.as.dataframe <- function(result.list){
	
	res = result.list$fitted
	
	lossP = sapply(res, "[[", "minLossPeriodic")
	lossNP = sapply(res, "[[", "minLossNonPeriodic")
	score = sapply(res, "[[", "score")
	ids = sapply(res, "[[", "ID")
	phi = sapply(res, "[[", "phi")
	lambda = sapply(res, "[[", "lambda")
	sigma = sapply(res, "[[", "sigma")
	amplitude = round(sapply(res, "[[", "a.coef"),2)
	mean = round(sapply(res, "[[", "b.coef"),2)
	psi = sapply(res, "[[", "psi")
	
	if(is.null(unlist(psi))){
		df = data.frame(ID=ids,score=round(score,2),phi=phi,lambda=lambda,sigma=sigma,mean=mean,amplitude=amplitude)
	} else {
		df = data.frame(ID=ids,score=round(score,2),phi=rep(NA,length(ids)),lambda=lambda,sigma=sigma,mean=mean,amplitude=amplitude)
		for(r in 1:(dim(psi)[1])){
			df = cbind(df,round(psi[r,],3))
			colnames(df)[dim(df)[2]] = paste("psi",r,sep="")
		}		
	}
	return(df)
}


