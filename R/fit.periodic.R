fit.periodic <- function(mat,
		timepoints=NULL,
		phi=NULL,
		lambda=NULL,
		sigma=NULL, 
		psi=NULL,  
		weights=NULL
){	
		
	if(class(mat) == "ExpressionSet") {mat = exprs(mat)}
	
	###################################################### 
	if(!is.null(weights)){
		if(!length(mat) == length(weights)){ stop("Error: wrong dimension of weights matrix")}
		if(is.null(rownames(mat)) | is.null(rownames(weights)) | !rownames(mat) == rownames(weights)){
			message("Warning: rownames of weight matrix not equal to rownames of data matrix.")
			message("replacing rownames of both matrices with ID_1,ID_2 ... ID_n")
			rownames(mat) = paste("ID",1:dim(mat)[1],sep="_")
			rownames(weights) = paste("ID",1:dim(weights)[1],sep="_")
		}
	}
	
	if(is.null(dim(mat))){mat = matrix(mat,nrow=1)}
	
	if(is.null(rownames(mat))){
		rownames(mat) = paste("ID",1:dim(mat)[1],sep="_")
	}
	if(is.null(timepoints)){
		timepoints = 1:ncol(mat)
	}
	
	if(is.null(phi)){
		phi = timepoints
		n.phigrid = length(phi)
	}
	
	if(is.null(lambda)){
		lambda = timepoints
		n.lambdagrid = length(lambda)
	}
	
	if(is.null(sigma)){
		sigma = 0
		n.sigmagrid = length(sigma)
	}
	
	if(!is.null(psi)){
		n.unitgrid = psi +1 
		unitgrid = seq(from=0,to=2*pi,length=n.unitgrid+1)
		psigrid = matrix(unitgrid,nrow=1)
		
		for(j in 1:(psi-1)){
			hilf = apply(psigrid,2,extend,unitgrid)
			hilf = unlist(hilf)
			psigrid = matrix(hilf,nrow=j+1)
		}
		n.psigrid = ncol(psigrid)
	}
	
	lower.standard = timepoints[1]
	if(lower.standard == 0) lower.standard = 1
	upper.standard= timepoints[length(timepoints)-1]
	
	## non periodic test functions
	testfunmat.nonPeriodic = createLinearTimeCourses(timepoints)
		
	if(is.null(psi)){
		## periodic test functions without shape
		Glist = list()
		for (s in sigma){
			for (l in lambda){
				phi.l = phi[phi <= l]
				for (p in phi.l){
					Glist[[length(Glist)+1]] = list(	
							phi     = p,
							psi     = NULL,
							lambda  = l,
							sigma = s)		
				}
			}
		}
	} else {
		## periodic test functions with shape
		Glist = list()
		for (s in sigma){
			for (l in lambda){
				phi.l = phi[phi <= l]
				for (p in phi.l){
					for (ps in 1:ncol(psigrid)){
					Glist[[length(Glist)+1]] = list(	
							phi     = p,
							psi		= psigrid[,ps],
							lambda  = l,
							sigma = s)			
					}
				}
			}
		}
	}
	
	############################################################
	## computation time estimation	
	nFunctions = length(Glist)
	tests = nFunctions*dim(mat)[1]
	vec = mat[1,]
	
	if(nFunctions > 1000 & tests > 10000){
		
		message("Estimating time required to do the analysis .... ")
		# estimate time to create 1000 functions and perform 100 linear regressions
		testfmatrix = list()
		if(is.null(psi)){
			t1 = system.time((for(x in 1:1000){		
				testfmatrix[[length(testfmatrix)+1]] = makewave(Glist[[x]],timepoints,lower=lower.standard,upper=upper.standard,fun=fun.lognormal)
			}))
			t = t1[[1]] * nFunctions/1000
		} else {
			t1 = system.time((for(x in 1:100){
				testfmatrix[[length(testfmatrix)+1]] = makewave(Glist[[x]],timepoints,lower=lower.standard,upper=upper.standard,fun=fun.lognormal.psi,psi=Glist[[x]]$psi)
			}))
			t = t1[[1]] * nFunctions/100
		}
				
		t2 = system.time((for(x in 1:10000){	
			lm(vec ~ testfmatrix[[1]])
		}))

	cost = (t + t2[[1]] * tests/10000)/60
	message(paste(tests,"linear regressions will be performed. This will take approx.",round(cost,1),"minutes."))
	if(cost > 30){message("See the MoPS vignette for help on speeding things up.")}
	} else {message(paste(tests,"linear regressions will be performed. This will take < 5 minutes."))}
	
	############################################################
	## Creation of test functions
	testfmatrix = list()
	message("Creating test functions ....")
	message("0 %")
	percent = 10	
	if(!is.null(psi)){
		for(i in 1:length(Glist)){
			p = round(i/length(Glist)*100,0)
			if(p >= percent) {message(paste(p,"%"));percent=percent+10}
			testfmatrix[[length(testfmatrix)+1]] = makewave(Glist[[i]],timepoints,lower=lower.standard,upper=upper.standard,fun=fun.lognormal.psi,psi=Glist[[i]]$psi)
		}
	} else {
		for(i in 1:length(Glist)){
			p = round(i/length(Glist)*100,0)
			if(p >= percent) {message(paste(p,"%"));percent=percent+10}
			testfmatrix[[length(testfmatrix)+1]] = makewave(Glist[[i]],timepoints,lower=lower.standard,upper=upper.standard,fun=fun.lognormal)
		}	
	}
	###########################################################################
	## fitting to data
	ids = rownames(mat)
	n = length(ids)
		
	resultList = list()
	
	message("Fitting of test functions to data ....")
	message("0 %")
	percent = 10
	
	for(i in 1:n){
		p = round(i/n*100,0)
		if(p >= percent) {message(paste(p,"%"));percent=percent+10}
		
		vec = mat[i,] 
		# periodic
		all.losses = rep(NA,length(testfmatrix))
		for(j in 1:length(testfmatrix)){
			if(is.null(weights)){
				linmod = lm(vec ~ testfmatrix[[j]])
			} else{
				w = weights[i,]
				linmod = lm(vec ~ testfmatrix[[j]],weights=w)
			}
			loss = sum(linmod$residuals^2)
			R = cor(vec, testfmatrix[[j]],use="pairwise.complete.obs") 
			if(R < 0 | is.na(R)){loss = Inf}
			all.losses[j] = loss
			
		}
		minPos.periodic =  which.min(all.losses)
		loss.periodic =  all.losses[minPos.periodic]
		# calculation of the linear model for the best fitting test function
		if(is.null(weights)){
			linmod = lm(vec ~ testfmatrix[[minPos.periodic]])
		} else{
			linmod = lm(vec ~ testfmatrix[[minPos.periodic]],weights=w)
		}
		a = linmod$coefficients[[2]]
		b = linmod$coefficients[[1]]
		
		
		# non periodic
		lossL.nonPeriodic = c()
		for(j in 1:dim(testfunmat.nonPeriodic)[2]){
			if(is.null(weights)){
				if(sum(is.na(vec)) > 0){ 
					linmod = lm(vec ~ testfunmat.nonPeriodic[,j])
				} else {
					linmod = lm(vec ~ testfunmat.nonPeriodic[,j])
				}
			} else {
				linmod = lm(vec ~ testfunmat.nonPeriodic[,j],weights=w)
			}
			loss = sum(linmod$residuals^2)
			lossL.nonPeriodic = c(lossL.nonPeriodic,loss)
		}
		minPos.nonPeriodic = which.min(lossL.nonPeriodic)
		loss.nonPeriodic = lossL.nonPeriodic[minPos.nonPeriodic]
		if(is.null(weights)){
			linmod = lm(vec ~ testfunmat.nonPeriodic[,minPos.nonPeriodic])
		} else{
			linmod = lm(vec ~ testfunmat.nonPeriodic[,minPos.nonPeriodic],weights=w)
		}
		loss.nonPeriodic = sum(linmod$residuals^2)
		is.weakly.periodic = loss.periodic <= loss.nonPeriodic
		score = log(loss.nonPeriodic/loss.periodic)
		resultList[[i]] = list(ID=ids[i],score = score,
				minLossPeriodic=loss.periodic,
				minLossNonPeriodic=loss.nonPeriodic,
				phi=Glist[[minPos.periodic]]$phi,
				psi=Glist[[minPos.periodic]]$psi,
				lambda=Glist[[minPos.periodic]]$lambda,
				sigma=Glist[[minPos.periodic]]$sigma,
				a.coef = a,
				b.coef = b
		)
	}
	result = list(fitted=resultList,time=timepoints,cols.mat=ncol(mat),phi=phi,lambda=lambda,sigma=sigma)
	return(result)
}



