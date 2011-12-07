#mapFilterProbe
setMethod("initialize","mapFilterProbe",function(.Object,filteredIndices,chromosome,position,strand)
{
	.Object@filteredIndices <- filteredIndices
	.Object@chromosome <- chromosome
	.Object@position <- position
	.Object@strand <- strand
	return(.Object)
})

#WfmFit
setMethod("initialize","Wfm",function(.Object,betaMAP,varbetaMAP,smoothPar,varEps,dataOrigSpace,dataWaveletSpace,design,phenoData,method,genome.info,n.levels,probePosition,wave.filt,Kj,prior,alpha,delta,two.sided,rescale,sigProbes,regions,GlocRegions,FDR,CI,F,varF,eff,varEff)
{
	.Object@betaMAP <- betaMAP
	.Object@varbetaMAP <- varbetaMAP
	.Object@smoothPar <- smoothPar
	.Object@varEps <- varEps
	.Object@dataOrigSpace <- dataOrigSpace
	.Object@dataWaveletSpace <- dataWaveletSpace
	.Object@design <- design
	.Object@phenoData <- phenoData
	.Object@method <- method
	.Object@genome.info <- genome.info
	.Object@n.levels <- n.levels
	.Object@probePosition <- probePosition
	.Object@wave.filt <- wave.filt
	.Object@Kj <- Kj
	.Object@prior <- prior
	.Object@alpha <- alpha
	.Object@delta <- delta
	.Object@two.sided <- two.sided
	.Object@rescale <- rescale
	.Object@sigProbes <- sigProbes
	.Object@regions <- regions
	.Object@GlocRegions <- GlocRegions
	.Object@FDR <- FDR
	.Object@CI <- CI
	.Object@F <- F
	.Object@varF <- varF
	.Object@eff <- eff
	.Object@varEff <- varEff
	return(.Object)
})


#genomeInfo
setMethod("initialize","genomeInfo",function(.Object,chromosome,strand,minPos,maxPos)
{
	.Object@chromosome <- chromosome
	.Object@strand <- strand
	.Object@minPos <- minPos
	.Object@maxPos <- maxPos
	return(.Object)	
})

#waveTilingFeatureSet
setMethod("initialize","waveTilingFeatureSet",function(.Object,tfs)
{
	.Object <- tfs
	class(.Object) <- "waveTilingFeatureSet"
	attr(class(.Object),"package") <- "waveTiling"
	return(.Object)
})
