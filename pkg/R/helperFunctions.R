
## functions for data input

cel2TilingFeatureSet <- function(dataPath,annotationPackage)
{
	expFiles <- paste(dataPath,list.celfiles(dataPath),sep="/")
	tfs <- read.celfiles(expFiles,pkgname=annotationPackage)
	return(tfs)
}


extractChrFromFasta <- function(char,fasPath)
{
	splt <- unlist(strsplit(char,fasPath))
	hlp <- splt[seq(2,length(splt),2)]
	splt2 <- unlist(strsplit(hlp,".fas"))
	return(splt2)
}

pm2mm <- function(z)
{
	z <- complement(subseq(z,start=13,end=13))
	return(z)
}

## functions for model fitting

normVec <- function(vec)
{
	return(vec/sqrt(sum(vec^2)))
}

wave.transform <- function(x,n.levels=log(length(x),2),filt="haar") 
{
	N <- length(x)
	J <- n.levels
	if (N/2^J != trunc(N/2^J)) stop("Sample size is not divisible by 2^J")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf
	storage.mode(h) <- "double"
	g <- dict$lpf
	storage.mode(g) <- "double"
	start <- 1
	stop <- N/2
	Y <- rep(0,N)
	for (j in 1:J)
	{
		W <- V <- numeric(N/2^j)
		out <- .C("dwt", as.double(x), as.integer(N/2^(j - 1)),L, h, g, W =  as.double(W), V = as.double(V), PACKAGE = "waveslim")[6:7]
		Y[start:stop] <- out$W
		x <- out$V
		start <- start+N/(2^j)
		stop <- stop+N/2^(j+1)
	}
	stop <- N
	Y[start:stop] <- x
	return(Y)
}

wave.backtransform <- function(D,n.levels,filt="haar")
{
	p <- length(D)
	if ((p%%2^n.levels)!=0) stop("Sample size is not divisible by 2^n.levels")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf
	storage.mode(h) <- "double"
	g <- dict$lpf
	storage.mode(g) <- "double"
	stop <- p
	start <- p-p/2^n.levels+1
	X <- D[start:stop]
	for (j in n.levels:1)
	{
		stop <- start-1
		start <- stop-p/2^j+1
		N <- length(X)
		XX <- numeric(2*p/2^j)
		X <- .C("idwt", as.double(D[start:stop]), as.double(X), as.integer(N),L, h, g, out = as.double(XX), PACKAGE = "waveslim")$out
	}
	return(X)
}

wave.backtransformK <- function(Cum,n.levels,order=1,filt="haar")
{
	p <- length(Cum)
	if ((p%%2^n.levels)!=0) stop("Sample size is not divisible by 2^n.levels")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf^order
	storage.mode(h) <- "double"
	g <- dict$lpf^order
	storage.mode(g) <- "double"
	stop <- p
	start <- p-p/2^n.levels+1
	X <- Cum[start:stop]
	for (j in n.levels:1)
	{
		stop <- start-1
		start <- stop-p/2^j+1
		N <- length(X)
		XX <- numeric(2*p/2^j)
		X <- .C("idwt", as.double(Cum[start:stop]), as.double(X), as.integer(N),L, h, g, out = as.double(XX), PACKAGE = "waveslim")$out
	}
	return(X)
}

MapMar <- function(D,X,vareps,J=0,eqsmooth=TRUE)
{
	N <- nrow(X)
	q <- ncol(X)
	K <- ncol(D)
	if (J==0)
	{
		ends <- K
	} else
	{
	ends <- cumsum(K/2^c(1:J,J))
	}
	B <- matrix(0,nrow=q,ncol=K)
	varB <- B
	phi <- B
	if (eqsmooth)
	{
		out <- .C("MAPMARGEQSMOOTH",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	} else
	{
		out <- .C("MAPMARG",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	}
	names(out) <- c("beta_MAP","varbeta_MAP","phi")
	dim(out[[1]]) <- c(K,q)
	out[[1]] <- t(out[[1]])
	dim(out[[2]]) <- c(K,q)
	out[[2]] <- t(out[[2]])
	dim(out[[3]]) <- c(K,q)
	out[[3]] <- t(out[[3]])
	return(out)
}

MapMarImp <- function(D,X,vareps,J=0,eqsmooth=TRUE)
{
	N <- nrow(X)
	q <- ncol(X)
	K <- ncol(D)
	if (J==0)
	{
		ends <- K
	} else
	{
	ends <- cumsum(K/2^c(1:J,J))
	}
	B <- matrix(0,nrow=q,ncol=K)
	varB <- B
	phi <- B
	if (eqsmooth)
	{
		out <- .C("MAPMARGIMPEQSMOOTH",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	} else
	{
		out <- .C("MAPMARGIMP",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	}
	names(out) <- c("beta_MAP","varbeta_MAP","phi")
	dim(out[[1]]) <- c(K,q)
	out[[1]] <- t(out[[1]])
	dim(out[[2]]) <- c(K,q)
	out[[2]] <- t(out[[2]])
	dim(out[[3]]) <- c(K,q)
	out[[3]] <- t(out[[3]])
	return(out)
}


WaveMarEstVarJ <- function(Y,X,n.levels,wave.filt,prior=c("normal","improper"),eqsmooth=TRUE,saveall=FALSE,D,var.eps,max.it,tol=1e-6,trace=FALSE)
{
	N <- nrow(D)
	K <- ncol(D)		
	Kj <- K/2^c(1:n.levels,n.levels)
	ends <- cumsum(Kj)
	starts <- c(1,ends[-(n.levels+1)]+1)
	varEps <- rep(mad(D[,starts[1]:ends[1]])^2,n.levels+1)
	if (var.eps=="mad")
	{
		max.it <- -1
		saveall <- TRUE
	}
	if (max.it<0)
	{
		for (i in 2:(n.levels+1))
		{
			varEps[i] <- ifelse(i<(n.levels+1),mad(D[,starts[i]:ends[i]])^2,varEps[i-1])
		}
	}
	crit1 <- sum(D^2)
	if (max.it<1)
	{
		if (prior=="normal")
		{
			WaveFit <- MapMar(D,X,varEps,n.levels,eqsmooth)
		}
		if (prior=="improper")
		{
			WaveFit <- MapMarImp(D,X,varEps,n.levels,eqsmooth)
		}
	} else
	{
		if (trace==TRUE) cat("Gauss Seidel algorithm...")
		for (j in 1:max.it)
		{
			crit0 <- crit1
			if (prior=="normal")
			{
				WaveFit <- MapMar(D,X,varEps,n.levels,eqsmooth)
			}
			if (prior=="improper")
			{
				WaveFit <- MapMarImp(D,X,varEps,n.levels,eqsmooth)
			}
			for (i in 1:(n.levels+1))
			{
				VinvD <- D[,starts[i]:ends[i]]-X%*%(t(X)%*%D[,starts[i]:ends[i]]/(diag(t(X)%*%X)+1/WaveFit$phi[,starts[i]:ends[i]]))
				varEps[i] <- sum(D[,starts[i]:ends[i]]*VinvD)/(ends[i]-starts[i]+1)/N
			}
			crit1 <- -1/2*sum(log(t(X)%*%X%*%WaveFit$phi+rep(1,ncol(X))))-1/2*(log(varEps)%*%Kj)-N*K*log(2*pi)/2-N*K/2
			if (trace==TRUE) cat("\n iteration", j,"of",max.it)
			if ((abs((crit0-crit1)/crit0))<tol)
			{
				if (trace==TRUE) cat("\n Gauss Seidel algorithm converged \n")
				break
			}
		}
		if ((abs(crit0-crit1)/crit0)>tol) cat("\n Warning: The maximum number of iterations reached without convergence\n")
	}
	if (!saveall)
	{
		WaveFit$beta_MAP <- matrix()
		WaveFit$varbeta_MAP <- matrix()
		WaveFit$phi <- matrix()
	}
	WaveFit$varEps <- varEps
	WaveFit$n.levels <- n.levels
	WaveFit$K <- K
	WaveFit$N <- N
	WaveFit$Kj <- Kj
	WaveFit$X <- X
	WaveFit$wave.filt <- wave.filt
	return(WaveFit)
}

## plot
# annoFIle needs the following columns: "chromosome", "strand", "feature", "ID", "start", "end"
makeNewAnnotationTrack <- function(annoFile,chromosome,strand,minBase,maxBase,feature="exon",dp=NULL)
{
	if (is.null(dp))
	{
		dp <- DisplayPars(exon="lightblue")
	}
	return(makeAnnotationTrack(region=annoFile[(annoFile$chromosome==chromosome)&(annoFile$strand==strand)&(annoFile$start<maxBase)&(annoFile$end>minBase)&(annoFile$feature==feature),c("start","end","feature","group","ID")],dp=dp))
}


makeNewAnnotationTextOverlay <- function(annoFile,chromosome,strand,minBase,maxBase,region,feature=c("gene","transposable_element_gene"),y=.5,dp=NULL)
{
	if (is.null(dp))
	{
		dp=DisplayPars(cex=1)
	}
	annohlp <- annoFile[(annoFile$chromosome==chromosome)&(annoFile$strand==strand)&(annoFile$start<maxBase)&(annoFile$end>minBase)&(is.element(annoFile$feature,feature)),c("ID","start","end")]
	return(makeTextOverlay(annohlp$ID,x=(annohlp$start+annohlp$end)/2,y=y,region=region,dp=dp))
}

makeNewTranscriptRectangleOverlay <- function(sigRegions,locations,start,end,region=NULL,dp=NULL)
{
	if (is.null(dp))
	{
		dp=DisplayPars(color="black")
	}
	detectedRegionsSelect <- matrix(sigRegions[sigRegions[,1]<end&sigRegions[,2]>start,],ncol=2)
	detectedRegionsSelect[detectedRegionsSelect[,1]<start,1] <- start
	detectedRegionsSelect[detectedRegionsSelect[,2]>end,2] <- end
	return(makeRectangleOverlay(start=locations[detectedRegionsSelect[,1]],end =locations[detectedRegionsSelect[,2]],region=region,dp=dp)) 
}

#### makeDesign
makeDesign<-function(design=c("time","circadian","group","factorial"),replics, noGroups, factor.levels=NULL) {
# 	if (method=="twoGroup")
# 	{
# 		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=2)
# 		Xorig[,1] <- 1
# 		Xorig[,2] <- rep(c(1,-1),replics) ## contr.helmert
# 	}
	if (design=="time") #time
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=noGroups)
		orderedFactor <- factor(1:noGroups,ordered=TRUE)
		desPoly <- lm(rnorm(noGroups)~orderedFactor,x=TRUE)$x
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desPoly[,2:noGroups],2,rep,replics)
		
	}
	else if (design=="circadian") # circadian
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=3)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=3)
		Xorig[,1] <- 1
		Xorig[,2] <- rep(sin(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
		Xorig[,3] <- rep(cos(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
	}
	else  if (design=="group") # group
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=noGroups)
		desHelmert <- contr.helmert(noGroups)
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desHelmert[,1:(noGroups-1)],2,rep,replics)
	}
	else if (design=="factorial")
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
	else {
		stop("No proper design specified!")
	}	
	out <- list()
	X.qr <- qr(Xorig)
	out$Xorthnorm <- qr.Q(X.qr)
	out$Xorig <- Xorig
	return (out);
}

#### makeContrasts
makeContrasts <- function(contrasts, nlevels) {
	if (contrasts=="compare") {
		q <- nlevels*(nlevels-1)/2
		contr <- matrix(0,nrow=q,ncol=nlevels)
		hlp1 <- rep(2:nlevels,1:(nlevels-1))
		hlp2 <- unlist(sapply(1:(nlevels-1),function(x) seq(1:x)))
		for (i in 1:nrow(contr))
		{
			contr[i,hlp1[i]] <- 1
			contr[i,hlp2[i]] <- -1
		}
	}

	return (contr);
}


##### plot
plotWfm <- function(fit,inf,annoFile,minPos,maxPos,trackFeature="exon",overlayFeature=c("gene","transposable_element_gene"),two.strand=TRUE,plotData=TRUE,plotMean=TRUE,tracks=0)
{
	if (missing(annoFile)) {stop("Annotation File is missing!!")}
	Gloc <- getProbePosition(fit)
	if (missing(minPos)) {minPos<-min(Gloc)}
	if (missing(maxPos)) {maxPos<-max(Gloc)}
	chromosome <- getChromosome(fit)
	strand <- getStrand(fit)
	selID <- (1:length(Gloc))[Gloc>minPos & Gloc<maxPos]
	sta <- min(selID)
	end <- max(selID)
	minBase <- Gloc[sta]
	maxBase <- Gloc[end]
	trackCount <- 1
	trackInfo <- list()
	overlayInfo <- list()
	if (plotData==TRUE)
	{
		Y <- getDataOrigSpace(fit)
		replcs <- as.numeric(table(getPhenoInfo(fit)$group))
		trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(t(Y[,sta:end])),probeStart=Gloc[sta:end],dp=DisplayPars(color=rep(1:length(replcs),replcs),ylim=range(Y[,sta:end]),pointSize=.3,pch=sequence(replcs)))
		names(trackInfo)[trackCount] <- "Data"
		trackCount <- trackCount + 1
	}
	if (two.strand==TRUE)
	{
		trackInfo[[trackCount]] <- makeNewAnnotationTrack(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand="forward",feature=trackFeature,dp=NULL)
		names(trackInfo)[trackCount] <- "F"
		overlayInfo[[trackCount]] <- makeNewAnnotationTextOverlay(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand="forward",region=c(trackCount,trackCount),feature=overlayFeature,y=0.5)
		trackCount <- trackCount + 1
		trackInfo[[trackCount]] <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)
		trackCount <- trackCount + 1
		trackInfo[[trackCount]] <- makeNewAnnotationTrack(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand="reverse",feature=trackFeature)
		names(trackInfo)[trackCount] <- "R"
		overlayInfo[[trackCount]] <- makeNewAnnotationTextOverlay(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand="reverse",region=c(trackCount,trackCount),feature=overlayFeature,y=0.5)
		trackCount <- trackCount + 1
	} else
	{
		trackInfo[[trackCount]] <- makeNewAnnotationTrack(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand=strand,feature=trackFeature)
		if (strand=="forward")
		{
			names(trackInfo)[trackCount] <- "F"
		}
		if (strand=="reverse")
		{
			names(trackInfo)[trackCount] <- "R"
		}
		overlayInfo[[trackCount]] <- makeNewAnnotationTextOverlay(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand=strand,region=c(trackCount,trackCount),feature=overlayFeature,y=0.5)
		trackCount <- trackCount + 1
		gAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)
		trackCount <- trackCount + 1
	}
	effects <- getEff(inf)
	regions <- getRegions(inf)
	noGroups <- length(table(getPhenoInfo(fit)$group))
	if (inherits(inf,"WfmInfMeans")) {
		plotMean <- FALSE
	}
	if (plotMean==TRUE)
	{
		trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[1,sta:end])+c(-0.4,0.4),pointSize=.3,pch=1,lwd=1,type="line"))
		names(trackInfo)[trackCount] <- "Mean"
		overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[1]]),end(regions[[1]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
		trackCount <- trackCount + 1
	}
#	if (getWfmFitMethod(object)=="twoGroup" | getWfmFitMethod(object)=="circadian")
	if ((inherits(fit,"WfmFitFactor") & noGroups==2) | inherits(fit,"WfmFitCircadian"))  
	{
		if (tracks==0)
		{
			trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[2,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[2,sta:end])[1]),max(1.3,range(effects[2,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
			#if (getWfmFitMethod(object)=="twoGroup")
			if (inherits(inf,"WfmInfCompare")) 
			{
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[2,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[2,sta:end])[1]),max(1.3,range(effects[2,sta:end])[2])+c(-0.2,0.2)),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- "FC"
			} else
			{
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(sqrt(effects[2,sta:end]^2+effects[3,sta:end]^2)),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(sqrt(effects[2,sta:end]^2+effects[3,sta:end]^2))[1]),max(1.2,range(sqrt(effects[2,sta:end]^2+effects[3,sta:end]^2))))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- "Ampl"
			}
			overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[2]]),end(regions[[2]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
		}
	}
#	if (getWfmFitMethod(object)=="compareGroupsTime" | getWfmFitMethod(object)=="compareGroupsFactor")
	else if (inherits(inf,"WfmInfCompare"))
	{
		if (length(tracks)==1 & tracks[1]==0)
		{
			
			effectsToPlot <- rep(0,noGroups)
			effectId <- c(1,3,6,10,15)
			effectNames <- c("2-1","3-2","4-3","5-4","6-5")
			effectNo <- 1
			while (effectNo<length(effectsToPlot))
			{
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[effectId[effectNo]+1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[effectId[effectNo]+1,sta:end])[1]),max(1.3,range(effects[effectId[effectNo]+1,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- effectNames[effectNo]
				overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[effectId[effectNo]+1]]),end(regions[[effectId[effectNo]+1]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
				effectNo <- effectNo + 1
				trackCount <- trackCount + 1
			}
			trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[effectId[length(effectsToPlot)-2]+2,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[effectId[length(effectsToPlot)-2]+2,sta:end])[1]),max(1.3,range(effects[effectId[length(effectsToPlot)-2]+2,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
			names(trackInfo)[trackCount] <- "Last-First"
			overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[effectId[length(effectsToPlot)-2]+2]]),end(regions[[effectId[length(effectsToPlot)-2]+2]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
			
		}
		if (length(tracks)>1 | tracks[1]!=0)
		{
			if (length(tracks)>8)
			{
				stop("Too many tracks specified")
			}
			groupNames <- sapply(1:noGroups,function(x) paste("Gr",x,sep=""))
			#groupNames <- sapply(8:13,function(x) paste("D",x,sep=""))
			firstId <- rep(2:noGroups,1:(noGroups-1))
			lastId <- sequence(1:(noGroups-1))
			effectNames <- sapply(1:(noGroups*(noGroups-1)/2),function(x) paste(groupNames[firstId[x]],"-",groupNames[lastId[x]],sep=""))
			for (i in tracks)
			{
				#trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i+1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[i+1,sta:end])[1]),max(1.3,range(effects[i+1,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i+1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[,sta:end])+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- effectNames[i]
				overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[i+1]]),end(regions[[i+1]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
				trackCount <- trackCount + 1
			}
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
		}
		
	}
#	if (getWfmFitMethod(object)=="meansByGroupTime" | getWfmFitMethod(object)=="meansByGroupFactor")
	else if (inherits(inf,"WfmInfMeans"))
	{
		if (length(tracks)==1 & tracks[1]==0)
		{
			effectNames <- sapply(1:noGroups,function(x) paste("Gr",x,sep=""))
			for (i in 1:noGroups)
			{
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[i,sta:end])[1]),max(1.3,range(effects[i,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- effectNames[i]
				overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[i]]),end(regions[[i]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
				trackCount <- trackCount + 1
			}
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
		}
		if (length(tracks)>1 | tracks[1]!=0)
		{
			if (length(tracks)>8)
			{
				stop("Too many tracks specified")
			}
			#effectNames <- sapply(1:noGroups,function(x) paste("Gr",x,sep=""))
			effectNames <- paste("day",8:13,sep="")
			for (i in tracks)
			{
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[i,sta:end])[1]),max(1.3,range(effects[i,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				#trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[,sta:end])+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- effectNames[i]
				overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[i]]),end(regions[[i]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
				trackCount <- trackCount + 1
			}
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
		}
	}
#	if (getWfmFitMethod(object)=="effectsTime")
	else if (inherits(fit,"WfmFitTime") & inherits(inf,"WfmInfEffects"))
	{
		fct <- 0.2390457
		if (tracks==0 | (length(tracks)==2 & tracks[1]==1 & tracks[2]==2) | tracks==1 | tracks==2)
		{
			if (tracks==0)
			{
				tracks <- c(1,2)
			}
			effectNames <- c("Linear","Quadratic")
			for (i in tracks)
			{
				#trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i+1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[i+1,sta:end])+c(-0.5,0.5),pointSize=.3,pch=1,lwd=1,type="line"))
				trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[i+1,sta:end]*fct),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[i+1,sta:end]*fct)+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
				names(trackInfo)[trackCount] <- effectNames[i]
				overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[i+1]]),end(regions[[i+1]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
				trackCount <- trackCount + 1
			}
			hlp <- which(!unlist(lapply(overlayInfo,is.null)))
			if (two.strand==TRUE)
			{
				hlpAnno <- hlp[1:2]
			} else
			{
				hlpAnno <- hlp[1]
			}
			hlpAnnoIndex <- unlist(lapply(hlpAnno,function(x) length(overlayInfo[[x]]@xpos)>0))
			hlpEff <- hlp[!(hlp%in%hlpAnno)]
			hlpEffIndex <- unlist(lapply(hlpEff,function(x) length(overlayInfo[[x]]@start)>0))
			hlpIndex <- c(hlpAnnoIndex,hlpEffIndex)
			overlayInfoCorr <- lapply(hlp,function(x) overlayInfo[[x]])
			gdPlot(trackInfo,minBase=minBase,maxBase=maxBase,overlay=overlayInfoCorr[hlpIndex])
		} else
		{
				stop("tracks argument has been misspecified")
		}
	}
	else if (inherits(inf,"WfmInfCustom")) {
		stop("Custom Inference plot not implemented!")
	}
	else {
		stop("Unknown Design")	  
	}
}