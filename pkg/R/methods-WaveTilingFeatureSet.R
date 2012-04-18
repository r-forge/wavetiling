# temporary (?)
# I can only use this after inserting strand information to the database when parsing the bpmap and Cel file in the (updated) package PDInfoBuilder
# setMethod("pmStrand","WaveTilingFeatureSet",function(object)
# {
# 	conn <- db(object)
# 	sql <- paste("SELECT fid, strand", "FROM pmfeature", "INNER JOIN chrom_dict", "USING(chrom)")
# 	tmp <- dbGetQuery(conn, sql)
# 	tmp <- tmp[order(tmp[["fid"]]),]
#         return(tmp[["strand"]])
# }
# )


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

setMethod("filterOverlap",signature("WaveTilingFeatureSet"),function(object,remap=TRUE,BSgenomeObject,chrId,strand=c("forward","reverse","both"),MM=FALSE)
{
	if (!inherits(object,"WaveTilingFeatureSet")) #class(object)!="WaveTilingFeatureSet")
	{
		stop("class of object is not WaveTilingFeatureSet.")
	}
	pmIndex <- oligo::pmindex(object)
	dataPMSeq <- pmSequence(object)
	if (remap)
	{
		setTrustBand <- min(unique(width(dataPMSeq)))
		if (setTrustBand < 15)
		{
			warning("There are probes smaller than 15 bp")
		}
		minTrustBand <- 15
		trBand <- max(setTrustBand,minTrustBand)
		if ((strand=="forward") | (strand=="both"))
		{
			message("begin remapping forward strand...\n")
			message("extract probe sequences\n")
			pmSeqDict <- PDict(dataPMSeq,tb.start=1,tb.end=trBand)
			message("match probe sequences to DNA sequence\n")
			chrSeqList <- list()
			for (i in chrId)
			{
				chrSeqList[[i]] <- BSgenomeObject[[i]]
			}
			names(chrSeqList) <- seqnames(BSgenomeObject)[chrId]
			startPM <- lapply(chrSeqList,function(x) startIndex(matchPDict(pmSeqDict,x)))
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
			message("remapping forward strand done\n")
		}
		if ((strand=="reverse") | (strand=="both"))
		{
			message("begin remapping reverse strand...\n")
			## Fix me: sometimes reverseComplement/sometimes not?
			dataPMSeqRevComp <- reverseComplement(pmSequence(object))
			#dataPMSeqRevComp <- pmSequence(object)
			cat("extract probe sequences\n")
			pmSeqDictRevComp <- PDict(dataPMSeqRevComp,tb.start=1,tb.end=trBand)
			if (strand=="reverse")
			{
				chrSeqList <- list()
				for (i in chrId)
				{
					chrSeqList[[i]] <- BSgenomeObject[[i]]
				}
				names(chrSeqList) <- seqnames(BSgenomeObject)[chrId]
			}
			message("match probe sequences to DNA sequence\n")
			startPMRevComp <- lapply(chrSeqList,function(x) startIndex(matchPDict(pmSeqDictRevComp,x)))
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
			message("remapping reverse strand done\n")
		}
	}
	else
	{
		if ((strand=="forward") | (strand=="both"))
		{
			pmMatchIndex <- which(pmStrand(object)==1 & pmChr(object) %in% seqnames(BSgenomeObject)[chrId])
			## Fix me: Is it always capital Chr?
			# sometimes 0/1, sometimes +/- ?
			chrInit <- pmChr(object)[pmMatchIndex]
			posInit <- pmPosition(object)[pmMatchIndex]
			strandInit <- pmStrand(object)[pmMatchIndex]
		}
		if ((strand=="reverse") | (strand=="both"))
		{
			dataPMSeqRevComp <- reverseComplement(pmSequence(object))
			pmMatchIndexRevComp <- which(pmStrand(object)==0 & pmChr(object) %in% seqnames(BSgenomeObject)[chrId])
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
		message("filter overlaps forward-reverse strand\n")
		## should be done more efficiently with something in Biostrings... ("%in%" does not work)
		dataPMSeqChar <- as.character(dataPMSeq)
		dataPMSeqRevCompChar <- as.character(dataPMSeqRevComp)
		forwardReverseOverlap <- which(dataPMSeqChar %in% dataPMSeqRevCompChar)
		reverseForwardOverlap <- which(dataPMSeqRevCompChar %in% dataPMSeqChar)
	}
	if (MM)
	{
		message("filter overlaps PM-MM\n")
		if (strand=="forward")
		{
			dataMMSeq <- pm2mm(dataPMSeq)
			dataPMSeqChar <- as.character(dataPMSeq)
			dataMMSeqChar <- as.character(dataMMSeq)
			pmMmOverlap <- which(dataPMSeqChar %in% dataMMSeqChar)
		}
		if (strand=="reverse")
		{
			dataMMSeqRevComp <- pm2mm(dataPMSeqRevComp)
			dataPMSeqRevCompChar <- as.character(dataPMSeqRevComp)
			dataMMSeqRevCompChar <- ac.character(dataMMSegRevComp)
			pmMmOverlap <- which(dataPMSeqRevCompChar %in% dataMMSeqRevCompChar)
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
	out <- new("MapFilterProbe",filteredIndices=filteredIndices,chromosome=filteredChromosome,position=filteredPosition,strand=filteredStrand)
	return(out)
}
)

setMethod("selectProbesFromTilingFeatureSet",signature("WaveTilingFeatureSet"),function(object,chromosome,strand=c("forward","reverse"),minPos,maxPos)
{
	if (strand=="forward")
	{
		strand <- 1
	}
	if (strand=="reverse")
	{
		strand <- 0
	}
	if (!inherits(object,"TilingFeatureSet")) #class(object)!="TilingFeatureSet")
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
	if (maxPos > max(pmPosition(object)))
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



setMethod("wfm.fit",signature("WaveTilingFeatureSet"),function(object,filter.overlap=NULL,design=c("time","circadian","group","factorial","custom"),n.levels,factor.levels=NULL,chromosome,strand,minPos,maxPos,design.matrix=NULL,var.eps=c("margLik","mad"),prior=c("normal","improper"),eqsmooth=FALSE,max.it=20,wave.filt="haar",skiplevels=NULL,trace=FALSE,save.obs=c("plot","regions","all"))
{
# construct filtered data set
	if ((names(pData(object))[1]!="group")|((names(pData(object))[2]!="replicate")))
	{
		stop("phenoData of WaveTilingFeatureSet object is not correct. Use 'addPheno()' first.")
	}
	if (!is.null(filter.overlap))
	{
		if (!inherits(filter.overlap,"MapFilterProbe"))
		{
			stop("class of filter.overlap is not MapFilterProbe. Use 'filterOverlap()' to create such an object.")
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
		dataInit <- data.frame(cbind(exprs(object)[probeId$selectionInit,],getPosition(filter.overlap)[probeId$selectionInit]))
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
		if (design=="custom") {
		    stop("Missing design.matrix!")
		}
		## here we make use of makeDesign and method
		desgn <- makeDesign(design=design,replics=getReplics(object),noGroups=noGroups,factor.levels=factor.levels)
		Z <- desgn$Xorthnorm
		X <- desgn$Xorig
	} else
	{
		if (mode(design.matrix) != "numeric") stop("design.matrix must be a numeric matrix")
		## must matrix be full rank?
		X <- design.matrix
		Z <- qr.Q(qr(X))
	}

	N <- nrow(Y)
	D <- t(apply(Y,1,wave.transform,n.levels,filt=wave.filt))
	fit <- WaveMarEstVarJ(Y=Y,X=Z,n.levels=n.levels,wave.filt=wave.filt,D=D,var.eps=var.eps,prior=prior,eqsmooth=eqsmooth,max.it=max.it,tol=1e-6,trace=trace,saveall=TRUE)	

	if (design=="circadian") {
	      noBetas<-3	
	}
	else {
	      noBetas <- noGroups
	}

	if (is.null(skiplevels))
	{
		skiplevels <- rep(0,noBetas)
	} else if (length(skiplevels)==1)
	{
		skiplevels <- rep(skiplevels,noBetas)
	}
	F <- matrix(0,nrow=noBetas,ncol=P)
 	varF <- matrix(0,nrow=noBetas,ncol=P)

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

## voorlopig laten staan
	if (save.obs=="plot")
	{
		fit$beta_MAP <- matrix()
		fit$varbeta_MAP <- matrix()
		fit$phi <- matrix()
	}

	genomeLoc <- new("GenomeInfo",chromosome=chromosome,strand=strand,minPos=min(Gloc),maxPos=max(Gloc))
        replics <- getReplics(object)

        if (design=="time") {
	  fitObject <- new("WfmFitTime",betaWav=fit$beta_MAP,varbetaWav=fit$varbeta_MAP,smoothPar=fit$phi,varEps=fit$varEps,dataOrigSpace=Y,dataWaveletSpace=D,design.matrix=X,phenoData=pData(object),genome.info=genomeLoc,n.levels=n.levels,probePosition=Gloc,wave.filt=wave.filt,Kj=fit$Kj,prior=prior,F=F,varF=varF,P=P,Z=Z,noGroups=noGroups,replics=replics)
	}
	else if (design=="circadian") {
	  fitObject <- new("WfmFitCircadian",betaWav=fit$beta_MAP,varbetaWav=fit$varbeta_MAP,smoothPar=fit$phi,varEps=fit$varEps,dataOrigSpace=Y,dataWaveletSpace=D,design.matrix=X,phenoData=pData(object),genome.info=genomeLoc,n.levels=n.levels,probePosition=Gloc,wave.filt=wave.filt,Kj=fit$Kj,prior=prior,F=F,varF=varF,P=P,Z=Z,noGroups=noGroups,replics=replics)
	}
	else if (design %in% c("group","factorial")) {
	  fitObject <- new("WfmFitFactor",betaWav=fit$beta_MAP,varbetaWav=fit$varbeta_MAP,smoothPar=fit$phi,varEps=fit$varEps,dataOrigSpace=Y,dataWaveletSpace=D,design.matrix=X,phenoData=pData(object),genome.info=genomeLoc,n.levels=n.levels,probePosition=Gloc,wave.filt=wave.filt,Kj=fit$Kj,prior=prior,F=F,varF=varF,P=P,Z=Z,noGroups=noGroups,replics=replics)
	}
	else if (design=="custom") {
	  fitObject <- new("WfmFitCustom",betaWav=fit$beta_MAP,varbetaWav=fit$varbeta_MAP,smoothPar=fit$phi,varEps=fit$varEps,dataOrigSpace=Y,dataWaveletSpace=D,design.matrix=X,phenoData=pData(object),genome.info=genomeLoc,n.levels=n.levels,probePosition=Gloc,wave.filt=wave.filt,Kj=fit$Kj,prior=prior,F=F,varF=varF,P=P,Z=Z,noGroups=noGroups,replics=replics)
	}

	return (fitObject);
})

