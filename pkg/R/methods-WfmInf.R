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
 	return(object@vareff)
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

setMethod("plot",signature=c(fit="WfmFit",inf="WfmInf"),plotWfm) 


