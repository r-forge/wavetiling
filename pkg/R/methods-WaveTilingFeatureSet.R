setMethod("addPheno",signature("WaveTilingFeatureSet"),function(object,noGroups,groupNames,replics,...)
{
	if (dim(object)[2] != sum(replics))
	{
		stop("Number of samples in object differs from total number of replicates.")
	}
	if (length(groupNames) != noGroups)
	{
		stop("Number of group names differs from number of groups")
	}
	if (length(replics)==1)
	{
		warning("Only one number of replicates given. This number is used for all groups")
		replics <- rep(replics,noGroups)
	}
	if ((length(replics)!=1) & (length(replics)!=noGroups))
	{
		stop("Length of replicates vector differs from number of groups")
	}
	pData <- data.frame(matrix(NA,sum(replics),2,dimnames=list(NULL,c("group","replicate"))))
	pData[1:sum(replics),1] <- rep(groupNames,replics)
	pData[1:sum(replics),2] <- sequence(replics)
	rownames(pData) <- paste(pData[,1],pData[,2],sep=".")
	metadata <- data.frame(labelDescription=c("sample group","replicate number within sample group"))
	phenoData <- new("AnnotatedDataFrame",data=pData,varMetadata=metadata)
	phenoData(object) <- phenoData
	return(object)
}
)

setMethod("getNoGroups",signature("WaveTilingFeatureSet"),function(object)
{
	if ((names(pData(object))[1]!="group")|((names(pData(object))[2]!="replicate")))
	{
		stop("phenoData of WaveTilingFeatureSet object is not correct. Use 'addPheno()' first.")
	}
	noGroups <- length(table(pData(object)$group))
	return(noGroups)
}
)

setMethod("getGroupNames",signature("WaveTilingFeatureSet"),function(object)
{
	if ((names(pData(object))[1]!="group")|((names(pData(object))[2]!="replicate")))
	{
		stop("phenoData of TilingFeatureSet object is not correct. Use 'addPheno()' first.")
	}
	groupNames <- unique(pData(object)$group)
	return(groupNames)
}
)

setMethod("getReplics",signature("WaveTilingFeatureSet"),function(object)
{
	if ((names(pData(object))[1]!="group")|((names(pData(object))[2]!="replicate")))
	{
		stop("phenoData of WaveTilingFeatureSet object is not correct. Use 'addPheno()' first.")
	}
	# keep the order safe
	orderHlp <- order(unique(pData(object)$group))
	replcsTemp <- as.numeric(table(pData(object)$group))
	replcs <- replcsTemp[orderHlp]
	return(replcs)
}
)

setMethod("filterOverlap",signature("WaveTilingFeatureSet"),function(object,remap=TRUE,fastaFile,chrId,strand=c("forward","reverse","both"),MM=FALSE)
{
	if (class(object)!="WaveTilingFeatureSet")
	{
		stop("class of object is not WaveTilingFeatureSet.")
	}
	pmIndex <- oligo::pmindex(object)
	dataPMSeq <- pmSequence(object)
	if (remap)
	{
		setTrustBand <- min(unique(nchar(dataPMSeq)))
		if (setTrustBand < 15)
		{
			warning("There are probes smaller than 15 bp")
		}
		minTrustBand <- 15
		trBand <- max(setTrustBand,minTrustBand)
		if ((strand=="forward") | (strand=="both"))
		{
			cat("begin remapping forward strand...\n")
			cat("extract probe sequences\n")
			pmSeqDict <- PDict(dataPMSeq,tb.start=1,tb.end=trBand)
			chrSeq <- readFASTA(fastaFile,strip.desc=FALSE)
			chrDNAString <- lapply(chrSeq,function(x) DNAString(paste(x[[2]],collapse="")))
			names(chrDNAString) <- paste("chr",chrId,sep="")
			cat("match probe sequences to DNA sequence\n")
			startPM <- lapply(chrDNAString,function(x) startIndex(matchPDict(pmSeqDict,x)))
			nposChrPM <- lapply(startPM,function(x) sapply(x,length))
			pmMatch <- Reduce("+",nposChrPM)==1
			pmMatchIndex <- (1:length(pmMatch))[pmMatch]
			index <- lapply(nposChrPM,function(x)
			{
				indexhlp <- x==1
				index <- indexhlp[pmMatchIndex]
			})
			mp <- stack(Map(which,index))
			chrInit <- c()
			chrInit[mp[,1]] <- as.character(mp[,2])
			posInit <- c()
			posInit[mp[,1]] <- unlist(lapply(startPM,function(x) unlist(x[pmMatch])))
			strandInit <- rep("forward",length(pmMatchIndex))
			cat("remapping reverse strand done\n")
		}
		if ((strand=="reverse") | (strand=="both"))
		{
			cat("begin remapping reverse strand...\n")
			## Fix me: sometimes reverseComplement/sometimes not?
			#dataPMSeqRevComp <- reverseComplement(pmSequence(object))
			dataPMSeqRevComp <- pmSequence(object)
			cat("extract probe sequences\n")
			pmSeqDictRevComp <- PDict(dataPMSeqRevComp,tb.start=1,tb.end=trBand)
			if (strand=="reverse")
			{
				
				chrSeq <- readFASTA(fastaFile,strip.desc=FALSE)
				chrDNAString <- lapply(chrSeq,function(x) DNAString(paste(x[[2]],collapse="")))
				names(chrDNAString) <- paste("chr",chrId,sep="")
			}
			cat("match probe sequences to DNA sequence\n")
			startPMRevComp <- lapply(chrDNAString,function(x) startIndex(matchPDict(pmSeqDictRevComp,x)))
			nposChrPMRevComp <- lapply(startPMRevComp,function(x) sapply(x,length))
			pmMatchRevComp <- Reduce("+",nposChrPMRevComp)==1
			pmMatchIndexRevComp <- (1:length(pmMatchRevComp))[pmMatchRevComp]
			index <- lapply(nposChrPMRevComp,function(x)
			{
				indexhlp <- x==1
				index <- indexhlp[pmMatchIndexRevComp]
			})
			mp <- stack(Map(which,index))
			chrInitRevComp <- c()
			chrInitRevComp[mp[,1]] <- as.character(mp[,2])
			posInitRevComp <- c()
			posInitRevComp[mp[,1]] <- unlist(lapply(startPMRevComp,function(x) unlist(x[pmMatchRevComp])))
			strandInitRevComp <- rep("reverse",length(pmMatchIndexRevComp))
			cat("remapping reverse strand done\n")
		}
	}
	else
	{
		if ((strand=="forward") | (strand=="both"))
		{
			pmMatch <- pmStrand(object)=="+"
			pmMatchIndex <- (1:length(pmMatch))[pmMatch]
			chrInit <- pmChr(object)[pmMatchIndex]
			posInit <- pmPosition(object)[pmMatchIndex]
			strandInit <- pmStrand(object)[pmMatchIndex]
		}
		if ((strand=="reverse") | (strand=="both"))
		{
			dataPMSeqRevComp <- reverseComplement(pmSequence(object))
			pmMatchRevComp <- pmStrand(object)=="-"
			pmMatchIndexRevComp <- (1:length(pmMatchRevComp))[pmMatchRevComp]
			chrInitRevComp <- pmChr(object)[pmMatchIndexRevComp]
			posInitRevComp <- pmPosition(object)[pmMatchIndexRevComp]
			strandInitRevComp <- pmStrand(object)[pmMatchIndexRevComp]
		}
	}
	forwardReverseOverlap <- NULL
	reverseForwardOverlap <- NULL
	pmMmOverlap <- NULL
	if (strand=="both")
	{
		cat("filter overlaps forward-reverse strand\n")
		forwardReverseOverlap <- (1:length(pmMatch))[(dataPMSeq %in% dataPMSeqRevComp)]
		reverseForwardOverlap <- (1:length(pmMatchRevComp))[(dataPMSeqRevComp %in% dataPMSeq)]
	}
	if (MM)
	{
		cat("filter overlaps PM-MM\n")
		if (strand=="forward")
		{
			dataMMSeq <- pm2mm(dataPMSeq)
			pmMmOverlap <- (1:length(pmMatch))[(dataPMSeq %in% dataMMSeq)]
		}
		if (strand=="reverse")
		{
			dataMMSeqRevComp <- pm2mm(dataPMSeqRevComp)
			pmMmOverlap <- (1:length(pmMatchRevComp))[(dataPMSeqRevComp %in% dataMMSeqRevComp)]
		}
	}
	indicesForward <- NULL
	indicesReverse <- NULL
	chromosomeForward <- NULL
	chromosomeReverse <- NULL
	positionForward <- NULL
	positionReverse <- NULL
	strandForward <- NULL
	strandReverse <- NULL
	if ((strand=="both") | (strand=="forward"))
	{
		indexClean <- (1:length(pmMatchIndex))[(!(pmMatchIndex %in% forwardReverseOverlap)) & (!(pmMatchIndex %in% pmMmOverlap))]
		indicesForward <- pmIndex[pmMatchIndex[indexClean]]
		chromosomeForward <- chrInit[indexClean]
		positionForward <- posInit[indexClean]
		strandForward <- strandInit[indexClean]
	}
	if ((strand=="both") | (strand=="reverse"))
	{
		indexRevCompClean <- (1:length(pmMatchIndexRevComp))[(!(pmMatchIndexRevComp %in% reverseForwardOverlap)) & (!(pmMatchIndexRevComp %in% pmMmOverlap))]
		indicesReverse <- pmIndex[pmMatchIndexRevComp[indexRevCompClean]]
		chromosomeReverse <- chrInitRevComp[indexRevCompClean]
		positionReverse <- posInitRevComp[indexRevCompClean]
		strandReverse <- strandInitRevComp[indexRevCompClean]
	}
	filteredIndices <- c(indicesForward,indicesReverse)
	filteredChromosome <- c(chromosomeForward,chromosomeReverse)
	filteredPosition <- c(positionForward,positionReverse)
	filteredStrand <- c(strandForward,strandReverse)
	filteredChromosome <- filteredChromosome[order(filteredIndices)]
	filteredPosition <- filteredPosition[order(filteredIndices)]
	filteredStrand <- filteredStrand[order(filteredIndices)]
	filteredIndices <- filteredIndices[order(filteredIndices)]
	out <- new("mapFilterProbe",filteredIndices=filteredIndices,chromosome=filteredChromosome,position=filteredPosition,strand=filteredStrand)
	return(out)
}
)

setMethod("selectProbesFromTilingFeatureSet",signature("WaveTilingFeatureSet"),function(object,chromosome,strand=c("forward","reverse"),minPos,maxPos)
{
	if (class(object)!="TilingFeatureSet")
		{
			stop("class of object is not TilingFeatureSet.")
		}
	if ((length(grep("chr",chromosome))>0) | (length(grep("Chr",chromosome))>0))
	{
		stop("give only the number (or letter) in the chromosome argument.")
	}
	if (missing(minPos))
	{
		minPos <- min(pmPosition(object))
	}
	if (minPos < min(pmPosition(object)))
	{
		minPos <- min(pmPosition(object))
	}
	if (missing(maxPos))
	{
		maxPos <- max(pmPosition(object))
	}
	if (maxPos < max(pmPosition(object)))
	{
		maxPos <- max(pmPosition(object))
	}
	if (minPos > maxPos)
	{
		stop("minPos is greater than maxPos")
	}
	selChrom <- (1:length(pmChr(object)))[pmChr(object)==paste("Chr",as.character(chromosome),sep="")]
	selStrand <- (1:length(pmStrand(object)))[pmStrand(object)==strand]
	# check what the outcome is of pmStrand() and adapt function accordingly
	selHlp <- intersect(selChrom,selStrand)
	selPos <- (1:length(pmPosition(object)))[(pmPosition(object)>=minPos)&(pmPosition(object)<=maxPos)]
	selection <- intersect(selHlp,selPos)
	return(selection)
}
)

setMethod("bgCorrQn",signature("WaveTilingFeatureSet"),function(object,useMapFilter=NULL)
{
	if (is.null(useMapFilter))
	{
		pmIndex <- pmindex(object)
	} else
	{
		pmIndex <- getFilteredIndices(useMapFilter)
	}
	sel <- exprs(object[pmIndex,])
	exprBgNorm <- normalize.quantiles(log(bg.adjust(sel),2))
	exprs(object) <- exprBgNorm
	return(object)
})


setMethod("makeDesign",signature("WaveTilingFeatureSet"),function(object,method=c("twoGroup","compareGroupsTime","compareGroupsFactor","circadian","twoFactors","meansByGroupTime","meansByGroupFactor","effectsTime"),factor.levels=NULL)
{
	replics <- getReplics(object)
	noGroups <- getNoGroups(object)
	if (method=="twoGroup")
	{
		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=2)
		Xorig[,1] <- 1
		Xorig[,2] <- rep(c(1,-1),replics)
	}
	if (method=="compareGroupsTime" | method=="meansByGroupTime" | method=="effectsTime")
	{
		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		orderedFactor <- factor(1:noGroups,ordered=TRUE)
		desPoly <- lm(rnorm(noGroups)~orderedFactor,x=TRUE)$x
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desPoly[,2:noGroups],2,rep,replics)
		
	}
	if (method=="circadian")
	{
		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=3)
		Xorig[,1] <- 1
		Xorig[,2] <- rep(sin(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
		Xorig[,3] <- rep(cos(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
	}
	if (method=="compareGroupsFactor" | method=="meansByGroupFactor")
	{
		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		desHelmert <- contr.helmert(noGroups)
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desHelmert[,1:(noGroups-1)],2,rep,replics)
	}
	if (method=="twoFactors")
	{
		if (is.null(factor.levels))
		{
			stop("No proper factor levels given.")
		}
		Xorig <- matrix(0,nrow=sum(replics),ncol=prod(factor.levels))
		Xorig[,1] <- 1
		desTrt1 <- contr.treatment(factor.levels[1])
		desTrt2 <- contr.treatment(factor.levels[2])
		#desHelmert1 <- contr.helmert(factor.levels[1])
		#desHelmert2 <- contr.helmert(factor.levels[2])
		replicsMat <- matrix(replics,nrow=factor.levels[1],ncol=factor.levels[2],byrow=TRUE)
		desFac1 <- apply(desTrt1,2,function(x) unlist(sapply(1:nrow(replicsMat),function(y) rep(x[y],sum(replicsMat[y,])))))
		Xorig[,2:(factor.levels[1])] <- desFac1
		desFac2 <- apply(desTrt2,2,function(x) unlist(sapply(1:nrow(replicsMat),function(y) rep(x,replicsMat[y,]))))
		Xorig[,(factor.levels[1]+1):(factor.levels[1]+factor.levels[2]-1)] <- desFac2
		# include interactions
		colno <- sum(factor.levels)
		for (i in 1:ncol(desFac1))
		{
			for (j in 1:ncol(desFac2))
			{
				Xorig[,colno] <- desFac1[,i]*desFac2[,j]
				colno <- colno + 1
			}
		}
	}
	out <- list()
	X.qr <- qr(Xorig)
	out$Xorthnorm <- qr.Q(X.qr)
	out$Xorig <- Xorig
	return(out)
})


setMethod("wfm.analysis",signature("WaveTilingFeatureSet"),function(object,filter.overlap=NULL,method=c("twoGroup","compareGroupsTime","compareGroupsFactor","circadian","meansByGroupTime","meansByGroupFactor","effectsTime","twoFactors"),n.levels,chromosome,strand,minPos,maxPos,design.matrix=NULL,var.eps=c("margLik","mad"),prior=c("normal","improper"),max.it=20,wave.filt="haar",skiplevels=NULL,trace=FALSE,save.obs=c("plot","regions","all"),alpha=0.05,nsim=1000,delta=NULL,rescale=NULL,two.sided=NULL,minRunPos=90,minRunProbe=1,factor.levels=NULL)
{
# construct filtered data set
	if ((names(pData(object))[1]!="group")|((names(pData(object))[2]!="replicate")))
	{
		stop("phenoData of WaveTilingFeatureSet object is not correct. Use 'addPheno()' first.")
	}
	
	if (!is.null(filter.overlap))
	{
		if (class(filter.overlap)!="mapFilterProbe")
		{
			stop("class of filter.overlap is not mapFilterProbe. Use 'filterOverlap()' to create such an object.")
		}
		if (missing(minPos))
		{
			minPos <- min(getPosition(filter.overlap))
		}
		if (missing(maxPos))
		{
			maxPos <- max(getPosition(filter.overlap))
		}
		probeId <- selectProbesFromFilterOverlap(filter.overlap,chromosome,strand,minPos=minPos,maxPos=maxPos)
		dataInit <- data.frame(cbind(exprs(object)[probeId$selectionFiltered,],getPosition(filter.overlap)[probeId$selectionFiltered]))
		#attention probeId$selection or probeId$selectionFiltered
		
	} else
	{
		if (missing(minPos))
		{
			minPos <- min(pmPosition(object))
		}
		if (missing(maxPos))
		{
			maxPos <- max(pmPosition(object))
		}
		probeId <- selectProbesFromTilingFeatureSet(object,chromosome,strand,minPos=minPos,maxPos=maxPos)
		dataInit <- data.frame(cbind(pm(object)[probeId,],pmPosition(object)[probeId]))
	}
	names(dataInit) <- c(rownames(pData(object)),"pos")
	# order expression values by genomic position
	dataInit <- dataInit[order(dataInit$pos),]
	if (missing(n.levels))
	{
		warning("Missing argument n.levels (levels of decomposition). Default value of 10 is used.")
		n.levels <- 10
	}
	nProbeInit <- dim(dataInit)[1]
	P <- nProbeInit - nProbeInit%%(2^n.levels)
	dataInit2 <- dataInit[1:P,]
	Y <- t(dataInit2[,1:(dim(dataInit2)[2]-1)])
	Gloc <- dataInit2[,"pos"]
	noGroups <- getNoGroups(object)
	# construct design matrix
	if (is.null(design.matrix))
	{
		desgn <- makeDesign(object=object,method=method,factor.levels=factor.levels)
		Z <- desgn$Xorthnorm
		X <- desgn$Xorig
	} else
	{
		X <- design.matrix
		Z <- qr.Q(qr(X))
	}
	N <- nrow(Y)
	D <- t(apply(Y,1,wave.transform,n.levels,filt=wave.filt))
	fit <- WaveMarEstVarJ(Y=Y,X=Z,n.levels=n.levels,wave.filt=wave.filt,D=D,var.eps=var.eps,prior=prior,max.it=max.it,tol=1e-6,trace=trace,saveall=TRUE)
	Xsel <- cumsum(getReplics(object))-getReplics(object)+1
	Xdes <- X[Xsel,]
	if (method=="twoGroup" | method=="compareGroupsTime" | method=="compareGroupsFactor")
	{
		q <- noGroups*(noGroups-1)/2
		noBetas <- noGroups
		contr <- matrix(0,nrow=q,ncol=noGroups)
		hlp1 <- rep(2:noGroups,1:(noGroups-1))
		hlp2 <- unlist(sapply(1:(noGroups-1),function(x) seq(1:x)))
		for (i in 1:nrow(contr))
		{
			contr[i,hlp1[i]] <- 1
			contr[i,hlp2[i]] <- -1
		}
		if (is.null(rescale))
		{
			rescale <- contr%*%Xdes
			rescale <- rbind(c(mean(Xdes[,1]),rep(0,noGroups-1)),rescale)
		}
	} else
	if (method=="effectsTime" | method=="twoFactors")
	{
		q <- noGroups-1
		noBetas <- noGroups
		rescale <- diag(noBetas)
	}
	if (method=="circadian")
	{
		noBetas <- 3
		q <- 1 # check
		rescale <- diag(noBetas)
	} else
	if (method=="meansByGroupTime" | method=="meansByGroupFactor")
	{
		q <- noGroups-1
		noBetas <- noGroups
		rescale <- Xdes
	}
	if (length(alpha)==1)
	{
		alpha <- rep(alpha,q+1)
	}
	F <- matrix(0,nrow=noBetas,ncol=P)
	varF <- matrix(0,nrow=noBetas,ncol=P)
	FDR <- matrix(0,nrow=q+1,ncol=P)
	CI <- rep(0,P*(q+1)*2)
	dim(CI) <- c(q+1,2,P)
	if (is.null(skiplevels))
	{
		skiplevels <- rep(0,noBetas)
	} else
	if (length(skiplevels)==1)
	{
		skiplevels <- rep(skiplevels,noBetas)
	}
	for (i in 1:noBetas)
	{
		if (skiplevels[i]>0)
		{
			F[i,] <- wave.backtransform(c(rep(0,sum(fit$Kj[1:skiplevels[i]])),fit$beta_MAP[i,-(1:(sum(fit$Kj[1:skiplevels[i]])))]),fit$n.levels ,filt=fit$wave.filt)
			varF[i,] <- wave.backtransformK(c(rep(0,sum(fit$Kj[1:skiplevels[i]])),fit$varbeta_MAP[i,-(1:(sum(fit$Kj[1:skiplevels[i]])))]),order=2,n.levels=fit$n.levels,filt=fit$wave.filt)
		} else
		{
			F[i,] <- wave.backtransform(fit$beta_MAP[i,],fit$n.levels ,filt=fit$wave.filt)
  			varF[i,] <- wave.backtransformK(fit$varbeta_MAP[i,],order=2,n.levels=fit$n.levels,filt=fit$wave.filt)
  		}
	}
	eff <- rescale%*%solve(t(Z)%*%X)%*%F
	varEff <- (rescale%*%solve(t(Z)%*%X))^2%*%varF
	sigProbes <- list()
	regions <- list()
	GlocRegions <- list()
	givenDelta <- delta
	if (is.null(givenDelta))
	{
		if (method=="meansByGroupTime" | method=="meansByGroupFactor")
		{
			delta <- rep(median(Y),q+1)
		} else
		{
			delta <- c(median(Y),rep(log(1.1,2),q))
		}
	} else
	if (length(givenDelta)==1)
	{
		delta <- rep(delta,q+1)
	} else
	if ((length(givenDelta)==2) & givenDelta[1]=="median")
	{
		delta <- rep(0,q+1)
		delta[1] <- median(Y)
		delta[2:(q+1)] <- rep(as.numeric(givenDelta[2]),q)
	} else
	if ((length(givenDelta)==q+1) & givenDelta[1]=="median")
	{
		delta <- rep(0,q+1)
		delta[1] <- median(Y)
		delta[2:(q+1)] <- as.numeric(givenDelta[2:(q+1)])
	}
	if (method=="twoGroup" | method=="compareGroupsTime" | method=="compareGroupsFactor" | method=="effectsTime" | method=="twoFactors")
	{
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
	} else
	if (method=="circadian")
	{
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
		nsim <- nsim
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
	} else
	if (method=="meansByGroupTime" | method=="meansByGroupFactor")
	{
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
	if (save.obs=="plot")
	{
		fit$beta_MAP <- matrix()
		fit$varbeta_MAP <- matrix()
		fit$phi <- matrix()
	}
	genomeLoc <- new("genomeInfo",chromosome=chromosome,strand=strand,minPos=min(Gloc),maxPos=max(Gloc))
	fitObject <- new("Wfm",betaMAP=fit$beta_MAP,varbetaMAP=fit$varbeta_MAP,smoothPar=fit$phi,varEps=fit$varEps,dataOrigSpace=Y,dataWaveletSpace=D,design=X,phenoData=pData(object),method=method,genome.info=genomeLoc,n.levels=n.levels,probePosition=Gloc,wave.filt=wave.filt,Kj=fit$Kj,prior=prior,alpha=alpha,delta=delta,two.sided=two.sided,rescale=rescale,sigProbes=sigProbes,regions=regions,GlocRegions=GlocRegions,F=F,varF=varF,eff=eff,varEff=varEff,FDR=FDR,CI=CI)
})
