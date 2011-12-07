setMethod("show",signature("mapFilterProbe"),function(object)
{
	cat("Remapped and filtered probe information\n")
	cat("No. of filtered probes:",length(object@filteredIndices),"\n")
}
)

setMethod("show",signature("Wfm"),function(object)
{
	cat("object of class 'Wfm'\n")
	cat("Fitted object from wavelet based functional model\n")
	cat("WFM method used:",object@method,"\n")
	cat("Wavelet filter used:",object@wave.filt,"\n")
	cat("Number of wavelet decomposition levels:",object@n.levels,"\n")
	cat("Number of probes used for estimation:",length(object@probePosition),"\n")
	cat("\n")
	cat("Genome info for fitted object:\n")
	cat("\tChromosome:",object@genome.info@chromosome,"\n")
	cat("\tStrand:",object@genome.info@strand,"\n")
	cat("\tMinimum probe position:",object@genome.info@minPos,"\n")
	cat("\tMaximum probe position:",object@genome.info@maxPos,"\n")
	cat("\n")
}
)
