setMethod("show",signature("MapFilterProbe"),function(object)
{
	cat("Remapped and filtered probe information\n")
	cat("No. of filtered probes:",length(object@filteredIndices),"\n")
}
)

setMethod("show",signature("GenomeInfo"),function(object)
{
	cat("Genome Info :\n")
	cat("\tChromosome:",object@chromosome,"\n")
	cat("\tStrand:",object@strand,"\n")
	cat("\tMinimum probe position:",format(object@minPos,scientific=FALSE),"\n")
	cat("\tMaximum probe position:",format(object@maxPos,scientific=FALSE),"\n")
	cat("\n")
}
)

setMethod("show",signature("WfmFit"),function(object)
{
	cat("Wavelet filter used:",object@wave.filt,"\n")
	cat("Number of wavelet decomposition levels:",object@n.levels,"\n")
	cat("Number of probes used for estimation:",length(object@probePosition),"\n")
	cat("\n")
	show(object@genome.info);
}
)

setMethod("show",signature("WfmFitTime"),function(object)
{
	cat("Fitted object from wavelet based functional model - Time Design\n")
	callNextMethod();
}
)

setMethod("show",signature("WfmFitFactor"),function(object)
{
	cat("Fitted object from wavelet based functional model - Factor/Group Design\n")
	callNextMethod();
}
)

setMethod("show",signature("WfmFitCircadian"),function(object)
{
	cat("Fitted object from wavelet based functional model - Circadian Design\n")
	callNextMethod();
}
)

setMethod("show",signature("WfmFitCustom"),function(object)
{
	cat("Fitted object from wavelet based functional model - Custom Design\n")
	callNextMethod();
}
)

setMethod("show",signature("WfmInf"),function(object)
{
	show(object@genome.info);
}
)

setMethod("show",signature("WfmInfCompare"),function(object)
{
	cat("Inference object from wavelet based functional model - Pairwise Comparison\n")
	callNextMethod();
})

setMethod("show",signature("WfmInfCustom"),function(object)
{
	cat("Inference object from wavelet based functional model - Custom Design\n")
	callNextMethod();
})

setMethod("show",signature("WfmInfEffects"),function(object)
{
	cat("Inference object from wavelet based functional model - Effects\n")
	cat("Time Effect if Time Design, Circadian Effect if Circadian Design\n")
	callNextMethod();
})

setMethod("show",signature("WfmInfMeans"),function(object)
{
	cat("Inference object from wavelet based functional model - Groupwise Means\n")
	callNextMethod();
})

setMethod("show",signature("WfmInfOverallMean"),function(object)
{
	cat("Inference object from wavelet based functional model - Overall Mean\n")
	callNextMethod();
})

