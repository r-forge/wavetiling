setMethod("getProbePosition",signature("WfmFit"),function(object)
{
	return(object@probePosition)	
})

setMethod("getNoProbes",signature("WfmFit"),function(object)
{
	return(length(object@probePosition))
})

setMethod("getBetaMAP",signature("WfmFit"),function(object)
{
	return(object@betaMAP)
})

setMethod("getVarBetaMAP",signature("WfmFit"),function(object)
{
	return(object@varbetaMAP)
})

setMethod("getSmoothPar",signature("WfmFit"),function(object)
{
	return(object@smoothPar)
})

setMethod("getVarEps",signature("WfmFit"),function(object)
{
	return(object@varEps)
})


setMethod("getGenomeInfo",signature("WfmFit"),function(object)
{
	return(object@genome.info)
})

setMethod("getChromosome",signature("WfmFit"),function(object)
{
	return(object@genome.info@chromosome)
})

setMethod("getStrand",signature("WfmFit"),function(object)
{
	return(object@genome.info@strand)
})

setMethod("getMinPos",signature("WfmFit"),function(object)
{
	return(object@genome.info@minPos)
})

setMethod("getMaxPos",signature("WfmFit"),function(object)
{
	return(object@genome.info@maxPos)
})

setMethod("getNoLevels",signature("WfmFit"),function(object)
{
	return(object@n.levels)
})


setMethod("getDesignMatrix",signature("WfmFit"),function(object)
{
	return(object@design.matrix)
})

setMethod("getPhenoInfo",signature("WfmFit"),function(object)
{
	return(object@phenoData)
})

setMethod("getDataOrigSpace",signature("WfmFit"),function(object)
{
	return(object@dataOrigSpace)
})

setMethod("getDataWaveletSpace",signature("WfmFit"),function(object)
{
	return(object@dataWaveletSpace)
})

setMethod("getWaveletFilter",signature("WfmFit"),function(object)
{
	return(object@wave.filt)
})

setMethod("getKj",signature("WfmFit"),function(object)
{
	return(object@Kj)
})

setMethod("getPrior",signature("WfmFit"),function(object)
{
	return(object@prior)
})


setMethod("getF",signature("WfmFit"),function(object)
{
	return(object@F)
})

setMethod("getP",signature("WfmFit"),function(object)
{
	return(object@P)
})

setMethod("getZ",signature("WfmFit"),function(object)
{
	return(object@Z)
})


setMethod("wfm.inference",signature("WfmFit"),function(object,contrast.matrix=NULL,contrasts=c("compare","means","effects","overallMean"),delta=NULL,two.sided=NULL,minRunPos=90,minRunProbe=1,alpha=0.05,nsim=1000,rescale=NULL)
{

	sigProbes <- list()
	regions <- list()
	GlocRegions <- list()
	givenDelta <- delta
	noGroups <- object@noGroups
	replics<- object@replics
	F<-object@F
	varF<-object@varF

	P<-ncol(F) ## so P does not need to be stored!

	Xsel <- cumsum(replics)-replics+1
	X <- object@design.matrix
	Xdes <- X[Xsel,]

	Z<-object@Z

	if (!is.null(contrast.matrix)) {
	      ## Given contrast matrix (CustomFit)
	      # Further implementation needed
	      warning("Custom Inference Procedure Not Implemented yet!")
	}
	else if (contrasts=="compare") {
		if (inherits(object,"WfmFitFactor") | inherits(object,"WfmFitTime") | inherits(object,"WfmFitCircadian" | inherits(object,"WfmFitCustom"))) {
			#q <- noGroups*(noGroups-1)/2
			q <- noGroups*(noGroups-1)/2
			contr <- makeContrasts(contrasts=contrasts,nlevels=noGroups);  
			noBetas <- noGroups
			if (is.null(rescale))
			{
				rescale <- contr%*%Xdes
				rescale <- rbind(c(mean(Xdes[,1]),rep(0,noGroups-1)),rescale)
			}
			## effects!
			eff <- rescale%*%solve(t(Z)%*%X)%*%F
			varEff <- (rescale%*%solve(t(Z)%*%X))^2%*%varF

			if (length(alpha)==1)
			{
				alpha <- rep(alpha,q+1)
			}
			FDR <- matrix(0,nrow=q+1,ncol=P)
			CI <- rep(0,P*(q+1)*2)
			dim(CI) <- c(q+1,2,P)
			if (is.null(givenDelta))
			{
				delta <- c(median(getDataOrigSpace(object)),rep(log(1.1,2),q))
			} else if (length(givenDelta)==1)
			{
				delta <- rep(delta,q+1)
			} else if ((length(givenDelta)==2) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- rep(as.numeric(givenDelta[2]),q)
			} else if ((length(givenDelta)==q+1) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- as.numeric(givenDelta[2:(q+1)])
			}

			if (is.null(two.sided))
			{
				two.sided <- c(0,rep(1,q))
			}
			for (i in 1:(q+1))
			{
				if (two.sided[i]==1)
				{
					#FDR[i,] <- pnorm(delta[i],abs(eff[i,]),sqrt(varEff[i,]))
					FDRUp <- pnorm(delta[i],eff[i,],sqrt(varEff[i,]))
					FDRDown <- 1-pnorm(-delta[i],eff[i,],sqrt(varEff[i,]))
					FDR[i,] <- pmin(FDRUp,FDRDown)
				}
				if (two.sided[i]==0)
				{
					FDR[i,] <- pnorm(delta[i],eff[i,],sqrt(varEff[i,]))
				}
				CI[i,1,] <- qnorm(alpha/2,eff[i,],sqrt(varEff[i,]))
				CI[i,2,] <- qnorm(1-alpha/2,eff[i,],sqrt(varEff[i,]))
			}
		} else {
			    stop ("Analysis 'compare' not specified for this design!")
		}
	} 
	else if (contrasts=="effects") {
		if (inherits(object,"WfmFitTime")) {
			#q <- noGroups-1
			q <- noGroups-1
			noBetas <- noGroups
			rescale <- diag(noBetas)
			if (length(alpha)==1)
			{
				alpha <- rep(alpha,q+1)
			}


		## effects!
			eff <- rescale%*%solve(t(Z)%*%X)%*%F
			varEff <- (rescale%*%solve(t(Z)%*%X))^2%*%varF
			FDR <- matrix(0,nrow=q+1,ncol=P)
			CI <- rep(0,P*(q+1)*2)
			dim(CI) <- c(q+1,2,P)
			if (is.null(givenDelta))
			{
				delta <- c(median(getDataOrigSpace(object)),rep(log(1.1,2),q))
			} else if (length(givenDelta)==1)
			{
				delta <- rep(delta,q+1)
			} else if ((length(givenDelta)==2) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- rep(as.numeric(givenDelta[2]),q)
			} else if ((length(givenDelta)==q+1) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- as.numeric(givenDelta[2:(q+1)])
			}

			if (is.null(two.sided))
			{
				two.sided <- c(0,rep(1,q))
			}
			for (i in 1:(q+1))
			{
				if (two.sided[i]==1)
				{
					#FDR[i,] <- pnorm(delta[i],abs(eff[i,]),sqrt(varEff[i,]))
					FDRUp <- pnorm(delta[i],eff[i,],sqrt(varEff[i,]))
					FDRDown <- 1-pnorm(-delta[i],eff[i,],sqrt(varEff[i,]))
					FDR[i,] <- pmin(FDRUp,FDRDown)
				}
				if (two.sided[i]==0)
				{
					FDR[i,] <- pnorm(delta[i],eff[i,],sqrt(varEff[i,]))
				}
				CI[i,1,] <- qnorm(alpha/2,eff[i,],sqrt(varEff[i,]))
				CI[i,2,] <- qnorm(1-alpha/2,eff[i,],sqrt(varEff[i,]))
			}
		}
		else if (inherits(object,"WfmFitCircadian")) {
			#q <- 1
			noBetas <- 3
			q <- 1 # check
			rescale <- diag(noBetas)
			if (length(alpha)==1)
			{
				alpha <- rep(alpha,q+1)
			}


		## effects!
			eff <- rescale%*%solve(t(Z)%*%X)%*%F
			varEff <- (rescale%*%solve(t(Z)%*%X))^2%*%varF
			FDR <- matrix(0,nrow=q+1,ncol=P)
			CI <- rep(0,P*(q+1)*2)
			dim(CI) <- c(q+1,2,P)
			if (is.null(givenDelta))
			{
				delta <- c(median(getDataOrigSpace(object)),rep(log(1.1,2),q))
			} else if (length(givenDelta)==1)
			{
				delta <- rep(delta,q+1)
			} else if ((length(givenDelta)==2) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- rep(as.numeric(givenDelta[2]),q)
			} else if ((length(givenDelta)==q+1) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- as.numeric(givenDelta[2:(q+1)])
			}

			## fix me: use of two.sided
			if (is.null(two.sided))
			{
				two.sided <- c(0,rep(1,q))
			}
			FDRUp <- pnorm(delta[1],eff[1,],sqrt(varEff[1,]))
			FDRDown <- 1-pnorm(-delta[1],eff[1,],sqrt(varEff[1,]))
			FDR[1,] <- pmin(FDRUp,FDRDown)
			CI[1,1,] <- qnorm(alpha/2,eff[1,],sqrt(varEff[1,]))
			CI[1,2,] <- qnorm(1-alpha/2,eff[1,],sqrt(varEff[1,]))
			FDR[2,] <- rep(0,P)
			
			for (k in 1:nsim)
			{
				effSinSampl <- rnorm(P,eff[2,],sqrt(varEff[2,]))
				effCosSampl <- rnorm(P,eff[3,],sqrt(varEff[3,]))
				amplSampl <- sqrt(effSinSampl^2 + effCosSampl^2)
				FDR[2,] <- FDR[2,] + (amplSampl < delta[2])
				if (k%%100==0) cat(k," ")
				# calculate CIs: to do
			}
			FDR[2,] <- FDR[2,]/nsim
		} 

		else {
			    stop ("Analysis 'effects' not specified for this design!")
		}
	} 
	 else if (contrasts=="means") {
		if (inherits(object,"WfmFitFactor") | inherits(object,"WfmFitTime") | inherits(object,"WfmFitCircadian") | inherits(object,"WfmFitCustom")) {
			#q <- noGroups-1
			q <- noGroups-1
			noBetas <- noGroups
			rescale <- Xdes
			if (length(alpha)==1)
			{
				alpha <- rep(alpha,q+1)
			}


		## effects!
			eff <- rescale%*%solve(t(Z)%*%X)%*%F
			varEff <- (rescale%*%solve(t(Z)%*%X))^2%*%varF
			FDR <- matrix(0,nrow=q+1,ncol=P)
			CI <- rep(0,P*(q+1)*2)
			dim(CI) <- c(q+1,2,P)
			if (is.null(givenDelta))
			{
				delta <- rep(median(getDataOrigSpace(object)),q+1)
			} else if (length(givenDelta)==1)
			{
				delta <- rep(delta,q+1)
			} else if ((length(givenDelta)==2) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- rep(as.numeric(givenDelta[2]),q)
			} else if ((length(givenDelta)==q+1) & givenDelta[1]=="median")
			{
				delta <- rep(0,q+1)
				delta[1] <- median(getDataOrigSpace(object))
				delta[2:(q+1)] <- as.numeric(givenDelta[2:(q+1)])
			}

			if (is.null(two.sided))
			{
				two.sided <- rep(1,q+1)
			}
			for (i in 1:(q+1))
			{
				FDR[i,] <- pnorm(delta[i],eff[i,],sqrt(varEff[i,]))
				CI[i,1,] <- qnorm(alpha/2,eff[i,],sqrt(varEff[i,]))
				CI[i,2,] <- qnorm(1-alpha/2,eff[i,],sqrt(varEff[i,]))
			}
		}
		else {
			    stop ("Analysis 'means' not specified for this design!")
		}
	}
	else if (contrasts == "overallMean") {
	      ## overallMean
	      # Further implementation needed
	      warning("Contrast 'overall mean' not yet implemented!")
	}
	else {
	      stop ("No contrast matrix of contrast statement specified!")
	}
	Gloc <- object@probePosition

	for (i in 1:(q+1))
	{
		sigProbes[[i]] <- (1:P)[FDR[i,]<alpha[i]]
		regions[[i]] <- cbind(sigProbes[[i]][c(TRUE,(sigProbes[[i]][2:length(sigProbes[[i]])]-sigProbes[[i]][1:(length(sigProbes[[i]])-1)])>1)],sigProbes[[i]][c((sigProbes[[i]][2:length(sigProbes[[i]])]-sigProbes[[i]][1:(length(sigProbes[[i]])-1)])>1,TRUE)])
		if (length(regions[[i]])==0)
		{
			regions[[i]] <- matrix(,0,2)
		}
		if (length(regions[[i]])==2)
		{
			regions[[i]] <- matrix(regions[[i]],nrow=1)
		}
		regions[[i]] <- matrix(regions[[i]][(Gloc[regions[[i]][,2]]-Gloc[regions[[i]][,1]])>minRunPos,],ncol=2)
		regions[[i]] <- matrix(regions[[i]][regions[[i]][,2]-regions[[i]][,1]>minRunProbe,],ncol=2)
		GlocRegions[[i]] <- matrix(cbind(Gloc[regions[[i]][,1]],Gloc[regions[[i]][,2]]),ncol=2)
		regions[[i]] <- IRanges(start=regions[[i]][,1],end=regions[[i]][,2])
		GlocRegions[[i]] <- IRanges(start=GlocRegions[[i]][,1],end=GlocRegions[[i]][,2])
	}
	#genomeLoc <- new("genomeInfo",chromosome=chromosome,strand=strand,minPos=min(Gloc),maxPos=max(Gloc))
	#infObject <- new("Wfm",method=method,alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI)
	if (!is.null(contrast.matrix))	{
		infObject <- new("WfmInfCustom",alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI,eff=eff,varEff=varEff)
	}
	else if (contrasts=="compare" ) {
		infObject <- new("WfmInfCompare",alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI,eff=eff,varEff=varEff)
	}
	else if (contrasts=="effects" ) {
		infObject <- new("WfmInfEffects",alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI,eff=eff,varEff=varEff)
	}
	else if (contrasts=="means" ) {
		infObject <- new("WfmInfMeans",alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI,eff=eff,varEff=varEff)
	}
	else if (contrasts=="overallMean" ) {
		infObject <- new("WfmInfOverallMean",alpha=alpha,delta=delta,two.sided=two.sided,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,FDR=FDR,CI=CI,eff=eff,varEff=varEff)
	}
	return (infObject);
})










