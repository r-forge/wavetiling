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

## remark: found this also as constructors for "NChannelSet","FeatureSet" and "TilingFeatureSet" (selectMethod("initialize","TilingFeatureSet"))
#waveTilingFeatureSet
setMethod("initialize","waveTilingFeatureSet",function (.Object, ...)
{
#    .local <- function (.Object, assayData, phenoData, ...)
#    {
#        mySlots <- slotNames(.Object)
#        dotArgs <- list(...)
#        isSlot <- names(dotArgs) %in% mySlots
#        if (missing(assayData)) {
#            assayData <- do.call(assayDataNew, dotArgs[!isSlot],
#                envir = parent.frame())
#        }
#        if (missing(phenoData)) {
#            phenoData <- annotatedDataFrameFrom(assayData, byrow = FALSE)
#        }
#        if (is.null(varMetadata(phenoData)[["channel"]])) {
#            varMetadata(phenoData)[["channel"]] <- factor(rep("_ALL_",
#                nrow(varMetadata(phenoData))), levels = c(assayDataElementNames(assayData),
#                "_ALL_"))
#        }
#        appl <- if (storageMode(assayData) == "list")
#            lapply
#        else eapply
#        assaySampleNames <- appl(assayData, function(elt) {
#            cnames <- colnames(elt)
#            if (is.null(cnames))
#                sampleNames(phenoData)
#            else cnames
#        })
#        sampleNames(assayData) <- assaySampleNames
#        sampleNames(phenoData) <- sampleNames(assayData)
#        do.call(callNextMethod, c(.Object, assayData = assayData,
#            phenoData = phenoData, dotArgs[isSlot]))
#    }
#    .local(.Object, ...)
	callNextMethod(.Object);
}
)

# > structure(function (.Object, ...)
# > {
# >     # check if the first argument is an ExpressionSet
# >     # if so: initialize the object with it and tells .local not to
# >     # overwrite slots corresponding to missing arguments.
# >     # otherwise: s
# >     overwrite.missing<- TRUE
# >     dotargs<- list(...)   
# >     if( length(dotargs)>  1&&  is(dotargs[[1]], 'ExpressionSet') ){       
# >         .Object<- dotargs[[1]]       
# >         overwrite.missing<- FALSE
# >         dotargs<- dotargs[-1]
# >     }       
# > 
# >     # .local should initialize (i.e. overwrite) a slot of .Object with
# > its prototype only if overwrite.missing=TRUE, or in any case with the
# > corresponding non missing argument for this slot.
# >      .local<- function (.Object, overwrite.missing=TRUE, assayData,
# > phenoData, featureData,
# >          exprs = new("matrix"), ...)
# >      {                  if (overwrite.missing&&  missing(assayData)) {
# >     # stuff ...
# >     }
# >     # other stuff ...
# >      }
# > 
# >     # call .local with overwrite.missing as its first argument
# >     do.call(.local, c(list(.Object, overwrite.missing), dotargs))
# >     
# > }

# setMethod("initialize","waveTilingFeatureSet",function(.Object, ...)
# {
# 	Object <- callNextMethod(.Object, ...)
# 	return(Object)
# })
