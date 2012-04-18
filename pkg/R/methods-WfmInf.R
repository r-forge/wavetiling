## accessors
 setMethod("getAlpha",signature("WfmInf"),function(object)
 {
 	return(object@alpha)
 })

 setMethod("getDelta",signature("WfmInf"),function(object)
 {
 	return(object@delta)
 })
 
 setMethod("getTwoSided",signature("WfmInf"),function(object)
 {
 	return(object@two.sided)
 })

 setMethod("getSigProbes",signature("WfmInf"),function(object)
 {
 	return(object@sigProbes)
 })
 
 setMethod("getRegions",signature("WfmInf"),function(object)
 {
 	return(object@regions)
 })
 
 setMethod("getGenomicRegions",signature("WfmInf"),function(object)
 {
 	return(object@GlocRegions)
 })
 
 setMethod("getFDR",signature("WfmInf"),function(object)
 {
 	return(object@FDR)
 })
 
 setMethod("getCI",signature("WfmInf"),function(object)
 {
 	return(object@CI)
 })

 setMethod("getEff",signature("WfmInf"),function(object)
 {
 	return(object@eff)
 })
 
 setMethod("getVarEff",signature("WfmInf"),function(object)
 {
 	return(object@varEff)
 })



setMethod("getSigGenes",signature(fit="WfmFit",inf="WfmInf"),function(fit,inf,annoFile)
{
	#Gloc <- getProbePosition(object)
	strand <- getStrand(fit)
	chromosome <- getChromosome(fit)
	regions <- getGenomicRegions(inf)
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
	message("find overlaps with detected regions...\n")
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


setMethod("getNonAnnotatedRegions",signature(fit="WfmFit",inf="WfmInf"),function(fit,inf,annoFile)
{
	#Gloc <- getProbePosition(object)
	strand <- getStrand(fit)
	chromosome <- getChromosome(fit)
	regions <- getGenomicRegions(inf)
	if (strand=="forward")
	{
		strandAlt <- "+"
		strandOpp <- "-"
	} else
	{
		strandAlt <- "-"
		strandOpp <- "+"
	}
	message("get annotated regions...\n")
	annoExons <- annoFile[(annoFile$strand==strandAlt)&(annoFile$chromosome==chromosome)&((annoFile$feature=="exon")|(annoFile$feature=="pseudogenic_exon")),c("chromosome","strand","feature","ID","start","end")]
	annoExonsOpp <- annoFile[(annoFile$strand==strandOpp)&(annoFile$chromosome==chromosome)&((annoFile$feature=="exon")|(annoFile$feature=="pseudogenic_exon")),c("chromosome","strand","feature","ID","start","end")]
	regGlocNoAnnoSense <- list()
	regGlocNoAnnoBoth <- list()
	nList <- length(regions)
	message("find overlaps with detected regions...\n")
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

setMethod("plotWfm",signature=c(fit="WfmFit",inf="WfmInf"),
function(fit,inf,annoFile,minPos,maxPos,trackFeature="exon",overlayFeature=c("gene","transposable_element_gene"),two.strand=TRUE,plotData=TRUE,plotMean=TRUE,tracks=0)
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
#       if (getWfmFitMethod(object)=="twoGroup" | getWfmFitMethod(object)=="circadian")
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
#       if (getWfmFitMethod(object)=="compareGroupsTime" | getWfmFitMethod(object)=="compareGroupsFactor")
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
#       if (getWfmFitMethod(object)=="meansByGroupTime" | getWfmFitMethod(object)=="meansByGroupFactor")
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
#       if (getWfmFitMethod(object)=="effectsTime")
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
)

