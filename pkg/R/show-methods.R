setMethod("show",signature("mapFilterProbe"),function(object)
{
	cat("Remapped and filtered probe information\n")
	cat("No. of filtered probes:",length(object@filteredIndices),"\n")
}
)

setMethod("show",signature("genomeInfo"),function(object)
{
	cat("Genome Info :\n")
	cat("\tChromosome:",object@chromosome,"\n")
	cat("\tStrand:",object@strand,"\n")
	cat("\tMinimum probe position:",object@minPos,"\n")
	cat("\tMaximum probe position:",object@maxPos,"\n")
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

# Needs further implementation :: what to show?
setMethod("show",signature("WfmInf"),function(object)
{
	cat("object of class 'WfmInf'\n")
	cat("\n")
}
)