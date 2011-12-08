setMethod("getProbePosition",signature("Wfm"),function(object)
{
	return(object@probePosition)	
})

setMethod("getNoProbes",signature("Wfm"),function(object)
{
	return(length(object@probePosition))
})

setMethod("getBetaMAP",signature("Wfm"),function(object)
{
	return(object@betaMAP)
})

setMethod("getVarBetaMAP",signature("Wfm"),function(object)
{
	return(object@varbetaMAP)
})

setMethod("getSmoothPar",signature("Wfm"),function(object)
{
	return(object@smoothPar)
})

setMethod("getVarEps",signature("Wfm"),function(object)
{
	return(object@varEps)
})


setMethod("getGenomeInfo",signature("Wfm"),function(object)
{
	return(object@genome.info)
})

setMethod("getChromosome",signature("Wfm"),function(object)
{
	return(object@genome.info@chromosome)
})

setMethod("getStrand",signature("Wfm"),function(object)
{
	return(object@genome.info@strand)
})

setMethod("getMinPos",signature("Wfm"),function(object)
{
	return(object@genome.info@minPos)
})

setMethod("getMaxPos",signature("Wfm"),function(object)
{
	return(object@genome.info@maxPos)
})

setMethod("getNoLevels",signature("Wfm"),function(object)
{
	return(object@n.levels)
})

setMethod("getWfmMethod",signature("Wfm"),function(object)
{
	return(object@method)
})

setMethod("getDesign",signature("Wfm"),function(object)
{
	return(object@design)
})

setMethod("getPhenoInfo",signature("Wfm"),function(object)
{
	return(object@phenoData)
})

setMethod("getDataOrigSpace",signature("Wfm"),function(object)
{
	return(object@dataOrigSpace)
})

setMethod("getDataWaveletSpace",signature("Wfm"),function(object)
{
	return(object@dataWaveletSpace)
})

setMethod("getWaveletFilter",signature("Wfm"),function(object)
{
	return(object@wave.filt)
})

setMethod("getKj",signature("Wfm"),function(object)
{
	return(object@Kj)
})

setMethod("getPrior",signature("Wfm"),function(object)
{
	return(object@prior)
})

setMethod("getAlpha",signature("Wfm"),function(object)
{
	return(object@alpha)
})

setMethod("getDelta",signature("Wfm"),function(object)
{
	return(object@delta)
})

setMethod("getTwoSided",signature("Wfm"),function(object)
{
	return(object@two.sided)
})

setMethod("getRescale",signature("Wfm"),function(object)
{
	return(object@rescale)
})

setMethod("getSigProbes",signature("Wfm"),function(object)
{
	return(object@sigProbes)
})

setMethod("getRegions",signature("Wfm"),function(object)
{
	return(object@regions)
})

setMethod("getGenomicRegions",signature("Wfm"),function(object)
{
	return(object@GlocRegions)
})

setMethod("getFDR",signature("Wfm"),function(object)
{
	return(object@FDR)
})

setMethod("getCI",signature("Wfm"),function(object)
{
	return(object@CI)
})

setMethod("getF",signature("Wfm"),function(object)
{
	return(object@F)
})

setMethod("getVarF",signature("Wfm"),function(object)
{
	return(object@varF)
})

setMethod("getEff",signature("Wfm"),function(object)
{
	return(object@eff)
})

setMethod("getVarEff",signature("Wfm"),function(object)
{
	return(object@vareff)
})

setMethod("plotWfm",signature("Wfm"),function(object,annoFile,minPos,maxPos,trackFeature="exon",overlayFeature=c("gene","transposable_element_gene"),two.strand=TRUE,plotData=TRUE,plotMean=TRUE,tracks=0)
# complete tracks option for different methods
{
	Gloc <- getProbePosition(object)
	chromosome <- getChromosome(object)
	strand <- getStrand(object)
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
		Y <- getDataOrigSpace(object)
		replcs <- as.numeric(table(getPhenoInfo(object)$group))
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
		overlayInfo[[trackCount]] <- 		 makeNewAnnotationTextOverlay(annoFile=annoFile,chromosome=chromosome,minBase=minBase,maxBase=maxBase,strand=strand,region=c(trackcount,trackCount),feature=overlayFeature,y=0.5)
		trackCount <- trackCount + 1
		gAxis <- makeGenomeAxis(add53 = TRUE,add35 = TRUE)
		trackCount <- trackCount + 1
	}
	effects <- getEff(object)
	regions <- getRegions(object)
	noGroups <- length(table(getPhenoInfo(object)$group))
	if (getWfmMethod(object)=="meansByGroupTime" | getWfmMethod(object)=="meansByGroupFactor")
	{
		plotMean <- FALSE
	}
	if (plotMean==TRUE)
	{
		trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[1,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=range(effects[1,sta:end])+c(-0.4,0.4),pointSize=.3,pch=1,lwd=1,type="line"))
		names(trackInfo)[trackCount] <- "Mean"
		overlayInfo[[trackCount]] <- makeNewTranscriptRectangleOverlay(sigRegions=as.matrix(data.frame(start(regions[[1]]),end(regions[[1]]))),location=Gloc,start=sta,end=end,region=c(trackCount,trackCount),dp=DisplayPars(color="darkgrey",alpha=.1))
		trackCount <- trackCount + 1
	}
	if (getWfmMethod(object)=="twoGroup" | getWfmMethod(object)=="circadian")
	{
		if (tracks==0)
		{
			trackInfo[[trackCount]] <- makeGenericArray(intensity=as.matrix(effects[2,sta:end]),probeStart=Gloc[sta:end],dp=DisplayPars(color="black",ylim=c(min(-1.3,range(effects[2,sta:end])[1]),max(1.3,range(effects[2,sta:end])[2]))+c(-0.2,0.2),pointSize=.3,pch=1,lwd=1,type="line"))
			if (getWfmMethod(object)=="twoGroup")
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
	if (getWfmMethod(object)=="compareGroupsTime" | getWfmMethod(object)=="compareGroupsFactor")
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
	if (getWfmMethod(object)=="meansByGroupTime" | getWfmMethod(object)=="meansByGroupFactor")
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
	if (getWfmMethod(object)=="effectsTime")
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
})


setMethod("getNonAnnotatedRegions",signature("Wfm"),function(object,annoFile)
{
	#Gloc <- getProbePosition(object)
	strand <- getStrand(object)
	chromosome <- getChromosome(object)
	regions <- getGenomicRegions(object)
	if (strand=="forward")
	{
		strandAlt <- "+"
		strandOpp <- "-"
	} else
	{
		strandAlt <- "-"
		strandOpp <- "+"
	}
	cat("get annotated regions...\n")
	annoExons <- annoFile[(annoFile$strand==strandAlt)&(annoFile$chromosome==chromosome)&((annoFile$feature=="exon")|(annoFile$feature=="pseudogenic_exon")),c("chromosome","strand","feature","ID","start","end")]
	annoExonsOpp <- annoFile[(annoFile$strand==strandOpp)&(annoFile$chromosome==chromosome)&((annoFile$feature=="exon")|(annoFile$feature=="pseudogenic_exon")),c("chromosome","strand","feature","ID","start","end")]
	regGlocNoAnnoSense <- list()
	regGlocNoAnnoBoth <- list()
	nList <- length(regions)
	cat("find overlaps with detected regions...\n")
	for (j in 1:nList)
	{
		regGlocIR <- regions[[j]]
		annoExonsIR <- IRanges(start=annoExons$start,end=annoExons$end)
		annoExonsOppIR <- IRanges(start=annoExonsOpp$start,end=annoExonsOpp$end)
		overSense <- findOverlaps(regGlocIR,annoExonsIR)
		overOpp <- findOverlaps(regGlocIR,annoExonsOppIR)
		noAnnoSenseId <- which(!(1:length(regGlocIR) %in% matchMatrix(overSense)[,1]))
		noAnnoOppId <- which(!(1:length(regGlocIR) %in% matchMatrix(overOpp)[,1]))
		noAnnoBothId <- which((1:length(regGlocIR) %in% noAnnoSenseId) & (1:length(regGlocIR) %in% noAnnoOppId))
		regGlocNoAnnoSense[[j]] <- regGlocIR[noAnnoSenseId]
		regGlocNoAnnoBoth[[j]] <- regGlocIR[noAnnoBothId]
	}
	noAnnoSenseIR <- regGlocNoAnnoSense[[1]]
	noAnnoBothIR <- regGlocNoAnnoBoth[[1]]
	for (i in 2:nList)
	{
		noAnnoSenseIRi <- regGlocNoAnnoSense[[i]]
		noAnnoSenseIR <- c(noAnnoSenseIR,noAnnoSenseIRi)
		noAnnoBothIRi <- regGlocNoAnnoBoth[[i]]
		noAnnoBothIR <- c(noAnnoBothIR,noAnnoBothIRi)
	}
	noAnnoSenseIRAll <- reduce(noAnnoSenseIR)
	noAnnoBothIRAll <- reduce(noAnnoBothIR)
	## TO DO:  include option to give maximum expression / FC per region
	out <- NULL
	out$noAnnoBoth <- GRanges(seqnames=Rle(rep(chromosome,length(noAnnoBothIRAll))),strand=Rle(rep(strandAlt,length(noAnnoBothIRAll))),ranges=noAnnoBothIRAll)
	out$noAnnoSense <- GRanges(seqnames=Rle(rep(chromosome,length(noAnnoSenseIRAll))),strand=Rle(rep(strandAlt,length(noAnnoSenseIRAll))),ranges=noAnnoSenseIRAll)
	return(out)
})


setMethod("getSigGenes",signature("Wfm"),function(object,annoFile)
{
	#Gloc <- getProbePosition(object)
	strand <- getStrand(object)
	chromosome <- getChromosome(object)
	regions <- getGenomicRegions(object)
	annoFile$strand[annoFile$strand=="forward"] <- "+"
	annoFile$strand[annoFile$strand=="reverse"] <- "-"
	annoFile$strand[!(annoFile$strand %in% c("+","-"))] <- "*"
	annoFileGR <- GRanges(seqnames=Rle(annoFile$chromosome),ranges=IRanges(start=annoFile$start,end=annoFile$end),strand=Rle(annoFile$strand),feature=annoFile$feature,id=annoFile$ID)
	if (strand=="forward")
	{
		strandAlt <- "+"
		strandOpp <- "-"
	} else
	{
		strandAlt <- "-"
		strandOpp <- "+"
	}
	annoChrGR <- annoFileGR[seqnames(annoFileGR)==chromosome]
	geneId <- which(values(annoChrGR)$feature=="gene" | values(annoChrGR)$feature=="transposable_element_gene")
	annoChrGeneGR <- annoChrGR[geneId]
	#annoChrGeneStrandGR <- annoChrGeneGR[strand(annoChrGeneGR)==strandAlt]
	#annoChrGeneStrandOppGR <- annoChrGeneGR[strand(annoChrGeneGR)==strandOpp]
	cat("find overlaps with detected regions...\n")
	nList <- length(regions)
	annoOver <- GRangesList()
	for (j in 1:nList)
	{
		regGlocIR <- regions[[j]]
		regGlocGR <- GRanges(seqnames=rep(chromosome,length(regGlocIR)),ranges=regGlocIR,strand=rep("*",length(regGlocIR)),effectNo=rep(j,length(regGlocIR)))
		overL <- findOverlaps(regGlocGR,annoChrGeneGR)
		regOver <- regGlocGR[queryHits(overL)]
		annoOverj <- annoChrGeneGR[subjectHits(overL)]
		overInt <- pintersect(regOver,annoOverj)
		values(annoOverj)$regNo <- queryHits(overL)
		values(annoOverj)$percOverGene <- width(overInt)/width(annoOverj)*100
		values(annoOverj)$percOverReg <- width(overInt)/width(regOver)*100
		totPercOverGeneHlp <- rep(0,max(subjectHits(overL)))
		totPercOverGeneHlp2 <- tapply(values(annoOverj)$percOverGene,subjectHits(overL),sum)
		totPercOverGeneHlp[as.numeric(names(totPercOverGeneHlp2))] <- totPercOverGeneHlp2
		values(annoOverj)$totPercOverGene <- totPercOverGeneHlp[subjectHits(overL)]
		annoOverj <- GRangesList(annoOverj)
		annoOver <- c(annoOver,annoOverj)
	}
	return(annoOver)
})

getATGPolynomWavelet <- function(bgCorrData,mapFilter,chromosome,strand,tair9GffFile,FCRegions,FCFit,thresh,aboveLim)
{
	thresh <- thresh
	threshdensity <- 0.005
	
	cat("Filter extreme outliers.\n")
	
	aboveid <- (1:dim(Y)[2])[apply(Y,2,function(x) max(x)>aboveLim)]
	above_IR <- IRanges(aboveid,aboveid)
	
	cat("Select differentially expressed genes (Compare).\n")
	
	FCpair_regions.mix <- FCRegions
	FCpair_EB_median <- FCFit
	FCpair_regionsCompare <- vector("list",15)
	spuriousID <- vector("list",15)
	filterID <- vector("list",15)
	ATG <- vector("list",15)
	for (j in 1:15)
	{
		FCpair_regions <- as.data.frame(FCpair_regions.mix[[j]])
		FCpairreg_probe <- IRanges(start=FCpair_regions[,1],end=FCpair_regions[,2])
		names(FCpair_regions) <- c("start","end")
		FCpair_regions$no_probes <- width(FCpairreg_probe)
		FCpair_regions[,1] <- Gloc[FCpair_regions.mix[[j]][,1]]
		FCpair_regions[,2] <- Gloc[FCpair_regions.mix[[j]][,2]]
		FCpair_regions$meanmedianFC <- sapply(1:dim(FCpair_regions)[1],function(x)
		{
			xx <- mean(FCpair_EB_median[j,FCpair_regions.mix[[j]][x,1]:FCpair_regions.mix[[j]][x,2]])
		})
		nreg <- dim(FCpair_regions)[1]
		FCpairreg <- IRanges(start=FCpair_regions[,1],end=FCpair_regions[,2])
		FCpair_regions$length <- width(FCpairreg)
		FCpair_regions$density <- FCpair_regions$no_probes/FCpair_regions$length
		exprOverlap_ana <- findOverlaps(t9chr_IR_expr_ana,FCpairreg)
		exprOverlap_ana150 <- findOverlaps(t9chr_IR_expr_ana,FCpairreg,maxgap=150)
		exprOverlap_opp <- findOverlaps(t9chr_IR_expr_opp,FCpairreg)
		exprnames <- t9chr_ana$attributes[exprID_ana[matchMatrix(exprOverlap_ana)[,1]]]
		exprnames_opp <- t9chr_opp$attributes[exprID_opp[matchMatrix(exprOverlap_opp)[,1]]]
		exon <- rep("",nreg)
		exonBis <- rep("",nreg)
		mRNA <- rep("",nreg)
		miRNA <- rep("",nreg)
		ncRNA <- rep("",nreg)
		pseudogenic_exon <- rep("",nreg)
		rRNA <- rep("",nreg)
		snoRNA <- rep("",nreg)
		snRNA <- rep("",nreg)
		transposable_element_gene <- rep("",nreg)
		tRNA <- rep("",nreg)
		feature <- rep("",nreg)
		opp_strand <- rep("",nreg)
		match150 <- rep(NA,nreg)
		noverlap <- length(exprnames)
		noverlap_opp <- length(exprnames_opp)
		for (i in 1:nreg)
		{
			match <- (1:noverlap)[matchMatrix(exprOverlap_ana)[,2]==i]
			matchFeature <- tair9gff_expr_ana$feature[matchMatrix(exprOverlap_ana)[match,1]]
			match_opp <- (1:noverlap_opp)[matchMatrix(exprOverlap_opp)[,2]==i]
			matchFeature_opp <- tair9gff_expr_opp$feature[matchMatrix(exprOverlap_opp)[match_opp,1]]
			exon[i] <- paste(unique(exprnames[match][matchFeature=="exon"]),collapse=", ")
			exonBis[i] <- paste(unique(exprnames[match][matchFeature=="exon"]),collapse="")
			opp_strand[i] <- paste(unique(exprnames_opp[match_opp][matchFeature_opp=="exon"]),collapse=", ")
			mRNA[i] <- paste(unique(exprnames[match][matchFeature=="mRNA"]),collapse=", ")
			miRNA[i] <- paste(unique(exprnames[match][matchFeature=="miRNA"]),collapse=", ")
			ncRNA[i] <- paste(unique(exprnames[match][matchFeature=="ncRNA"]),collapse=", ")
			pseudogenic_exon[i] <- paste(unique(exprnames[match][matchFeature=="pseudogenic_exon"]),collapse=", ")
			rRNA[i] <- paste(unique(exprnames[match][matchFeature=="rRNA"]),collapse=", ")
			snoRNA[i] <- paste(unique(exprnames[match][matchFeature=="snoRNA"]),collapse=", ")
			snRNA[i] <- paste(unique(exprnames[match][matchFeature=="snRNA"]),collapse=", ")
			transposable_element_gene[i] <- paste(unique(exprnames[match][matchFeature=="transposable_element_gene"]),collapse=", ")
			tRNA[i] <- paste(unique(exprnames[match][matchFeature=="tRNA"]),collapse=", ")
			featureID <- (1:9)[c(mRNA[i],miRNA[i],ncRNA[i],pseudogenic_exon[i],rRNA[i],snoRNA[i],snRNA[i],transposable_element_gene[i],tRNA[i])!=""]
			feature[i] <- ifelse(length(featureID)==0,"no strand annotation",paste(c("mRNA","miRNA","ncRNA","pseudogenic_exon","rRNA","snoRNA","snRNA","transposable_element_gene","tRNA")[featureID],collapse=", "))
			match150[i] <- i %in% matchMatrix(exprOverlap_ana150)[,2]
		}
		FCpair_regions <- data.frame(FCpair_regions,as.character(exon),as.character(feature),as.character(opp_strand),match150)
		names(FCpair_regions)[7:10] <- c("exon","feature","opp_strand","flankingAnno")
		spurious <- (1:nreg)[abs(FCpair_regions$meanmedian)<thresh]
		sparse <- (1:nreg)[FCpair_regions$density<threshdensity]
		overlapAbove <- findOverlaps(FCpairreg_probe,above_IR)
		Regoverlap <- FCpairreg_probe[matchMatrix(overlapAbove)[,1]]
		outid1 <- (1:length(Regoverlap))[width(Regoverlap)<9]
		outid2 <- matchMatrix(overlapAbove)[outid1,1]
		spsp <- unique(c(spurious,sparse,outid2))
		spuriousID[[j]] <- spsp
		#further prosprocessing
		FCpair_regionsNew <- FCpair_regions[-spsp,]
		FCpair_regionsNewIR <- IRanges(start=FCpair_regionsNew[,1],end=FCpair_regionsNew[,2])
		nReg <- length(FCpair_regionsNewIR)
		#select short regions
		lowwidth_id <- (1:nReg)[width(FCpair_regionsNewIR)<300]
		#select stand-alone regions
		hlp <- matchMatrix(findOverlaps(FCpair_regionsNewIR,FCpair_regionsNewIR,maxgap=300))
		nOver <- dim(hlp)[1]
		hlpid <- (1:nOver)[hlp[,1]!=hlp[,2]]
		notalone_id <- hlp[hlpid,1]
		alone_id <- (1:nOver)[!(1:nOver %in% notalone_id)]
		suspect_id <- (1:nOver)[(1:nOver %in% lowwidth_id)&(1:nOver %in% alone_id)]
		suspectReg <- FCpair_regionsNewIR[suspect_id]
		suspectProbe <- GReg2id(suspectReg,Gloc)
		gap <- rep(0,length(suspectProbe))
		for (i in 1:length(suspectProbe))
		{
			dat <- as.numeric(Y[,(max(1,(start(suspectProbe)[i]-8))):(min(end(suspectProbe)[i]+8,length(Gloc)))])
			ll <- length(dat)
			sdata <- sort(dat,decreasing=TRUE)
			gap[i] <- max(sdata[1:(ll-1)]-sdata[2:ll])
		}
		filter_id <- suspect_id[(1:length(suspectProbe))[gap > 1.5]]
		filterID[[j]] <- filter_id
		FCpair_regionsClean <- FCpair_regionsNew[-filter_id,]
		highFC_id <- order(abs(FCpair_regionsClean$meanmedianFC),decreasing=TRUE)
		FCpair_regionsOrd <- FCpair_regionsClean[highFC_id,]
		FCpair_regionsCompare[[j]] <- FCpair_regionsOrd
		exonBis <- exonBis[-spsp]
		exonBis <- exonBis[-filter_id]
		tt <- strsplit(exonBis,"Parent=")
		tt <- unlist(tt)
		tt <- tt[tt!=""]
		ATG[[j]] <- unique(tt)
	}
	for (i in 1:15)
	{
		FCpair_regions.mix[[i]] <- FCpair_regions.mix[[i]][-spuriousID[[i]],]
		FCpair_regions.mix[[i]] <- FCpair_regions.mix[[i]][-filterID[[i]],]
	}
	out$AtgCompareSplice <- ATG
	out$AdaptedFCPairRegions <- FCpair_regions.mix
	out$FCPairRegionsInfo <- FCpair_regionsCompare
	out$Gloc <- Gloc
	return(out)
}






