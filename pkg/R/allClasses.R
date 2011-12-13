#mapFilterProbe
setClass("mapFilterProbe",representation(filteredIndices="vector",chromosome="vector",position="vector",strand="vector"))

#genomeInfo
setClass("genomeInfo",representation(chromosome="vector",strand="character",minPos="numeric",maxPos="numeric"))

#Wfm
setClass("Wfm",representation(betaMAP="matrix",varbetaMAP="matrix",smoothPar="matrix",varEps="numeric",dataOrigSpace="matrix",dataWaveletSpace="matrix",design="matrix",phenoData="data.frame",method="character",genome.info="genomeInfo",n.levels="numeric",probePosition="vector",wave.filt="character",Kj="numeric",prior="character",alpha="numeric",delta="numeric",two.sided="numeric",rescale="matrix",sigProbes="list",regions="list",GlocRegions="list",FDR="matrix",CI="array",F="matrix",varF="matrix",eff="matrix",varEff="matrix"))

#WaveTilingFeatureSet
setClass("WaveTilingFeatureSet",contains="TilingFeatureSet")

