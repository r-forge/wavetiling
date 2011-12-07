setMethod("getFilteredIndices",signature("mapFilterProbe"),function(object)
{
	return(object@filteredIndices)
}
)

setMethod("getChromosome",signature("mapFilterProbe"),function(object)
{
	return(object@chromosome)
}
)

setMethod("getPosition",signature("mapFilterProbe"),function(object)
{
	return(object@position)
}
)

setMethod("getStrand",signature("mapFilterProbe"),function(object)
{
	return(object@strand)
}
)

setMethod("selectProbesFromFilterOverlap",signature("mapFilterProbe"),function(object,chromosome,strand=c("forward","reverse"),minPos=min(getPosition(object)),maxPos=max(getPosition(object)))
{
	if (class(object)!="mapFilterProbe")
		{
			stop("class of object is not mapFilterProbe. Use 'filterOverlap()' to create such an object.")
		}
	if ((length(grep("chr",chromosome))>0) | (length(grep("Chr",chromosome))>0))
	{
		stop("give only the number (or letter) in the chromosome argument.")
	}
	if (minPos > maxPos)
	{
		stop("minPos is greater than maxPos")
	}
	selChrom <- (1:length(getFilteredIndices(object)))[getChromosome(object)==paste("chr",as.character(chromosome),sep="")]
	selStrand <- (1:length(getFilteredIndices(object)))[getStrand(object)==strand]
	selHlp <- intersect(selChrom,selStrand)
	selPos <- (1:length(getFilteredIndices(object)))[(getPosition(object)>=minPos)&(getPosition(object)<=maxPos)]
	selectionInit <- intersect(selHlp,selPos)
	selection <- getFilteredIndices(object)[selectionInit]
	return(list(selection=selection,selectionFiltered=selectionInit))
}
)

