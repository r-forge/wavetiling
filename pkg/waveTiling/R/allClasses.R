#MapFilterProbe
setClass("MapFilterProbe",representation(filteredIndices="vector",chromosome="vector",position="vector",strand="vector"))

#GenomeInfo
setClass("GenomeInfo",representation(chromosome="vector",strand="character",minPos="numeric",maxPos="numeric"))

#WaveTilingFeatureSet
setClass("WaveTilingFeatureSet",contains="TilingFeatureSet")

#WfmFit
setClass(Class="WfmFit",
	representation = representation	(betaWav="matrix",varbetaWav="matrix",smoothPar="matrix",varEps="numeric",dataOrigSpace="matrix",dataWaveletSpace="matrix",design.matrix="matrix",phenoData="data.frame",genome.info="GenomeInfo",n.levels="numeric",probePosition="vector",wave.filt="character",Kj="numeric",prior="character",F="matrix",varF="matrix",P="numeric",Z="matrix",noGroups="numeric",replics="numeric")
)

### Factor
setClass(Class="WfmFitFactor", contains="WfmFit")

### Time
setClass(Class="WfmFitTime", contains="WfmFit")

## Circadian
setClass(Class="WfmFitCircadian", contains="WfmFit")

## Custom
setClass(Class="WfmFitCustom", contains="WfmFit")

#WfmInf
setClass(Class="WfmInf",
	representation = representation	(alpha="numeric",delta="numeric",two.sided="numeric",sigProbes="list",regions="list",GlocRegions="list",FDR="matrix",CI="array",eff="matrix",varEff="matrix")
)

## Compare
setClass(Class="WfmInfCompare", contains="WfmInf")

## Effects
setClass(Class="WfmInfEffects", contains="WfmInf")

## Means
setClass(Class="WfmInfMeans", contains="WfmInf")

## Overall Mean
setClass(Class="WfmInfOverallMean", contains="WfmInf")

## Custom
setClass(Class="WfmInfCustom", contains="WfmInf")
