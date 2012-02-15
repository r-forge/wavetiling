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
setMethod("initialize","WfmFit",function(.Object,betaMAP,varbetaMAP,smoothPar,varEps,dataOrigSpace,dataWaveletSpace,design.matrix,phenoData,genome.info,n.levels,probePosition,wave.filt,Kj,prior,F,varF,P,Z,noGroups,replics)
{
	.Object@betaMAP <- betaMAP
	.Object@varbetaMAP <- varbetaMAP
	.Object@smoothPar <- smoothPar
	.Object@varEps <- varEps
	.Object@dataOrigSpace <- dataOrigSpace
	.Object@dataWaveletSpace <- dataWaveletSpace
	.Object@design.matrix <- design.matrix
	.Object@phenoData <- phenoData
#	.Object@design <- design
	.Object@genome.info <- genome.info
	.Object@n.levels <- n.levels
	.Object@probePosition <- probePosition
	.Object@wave.filt <- wave.filt
	.Object@Kj <- Kj
	.Object@prior <- prior
#	.Object@rescale <- rescale
	.Object@F <- F
	.Object@varF <- varF
	.Object@P <- P
	.Object@Z <- Z
	.Object@noGroups <- noGroups
	.Object@replics <- replics
	return(.Object)
})


# WfmInf
setMethod("initialize","WfmInf",function(.Object, alpha, delta, two.sided, sigProbes, regions, GlocRegions, FDR, CI, eff, varEff)
{
	.Object@alpha <- alpha
	.Object@delta <- delta
	.Object@two.sided <- two.sided
	.Object@sigProbes <- sigProbes
	.Object@regions <- regions
	.Object@GlocRegions <- GlocRegions
	.Object@FDR <- FDR
	.Object@CI <- CI
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
setMethod("initialize","WaveTilingFeatureSet",function (.Object)
{
	callNextMethod(.Object);
}
)

## Subclasses
setMethod("initialize","WfmFitTime",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmFitCircadian",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmFitFactor",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmFitCustom",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmInfCustom",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmInfCompare",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmInfEffects",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmInfMeans",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)

setMethod("initialize","WfmInfOverallMean",function (.Object, ...)
{
	callNextMethod(.Object, ...);
}
)
