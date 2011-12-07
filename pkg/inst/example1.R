## example - read.celfiles ##
#############################

data_path = "../data/CEL"
#data_path <- "/home/kdebeuf/research/tiling/arabidopsis_data"
exp_names = list.celfiles(data_path)
dna.e2f.wt=read.celfiles(paste(data_path,exp_names,sep="/"),pkgname="pd.at35b.mr.v04.2.tigrv5")

#dna.e2f.wt is object of class TilingFeatureSet
class(dna.e2f.wt)
#[1] "TilingFeatureSet"
#attr(,"package")
#[1] "oligoClasses"

##Accessing Data Elements
#names of the features
featureNames(dna.e2f.wt)
#unique sample identifiers
sampleNames(dna.e2f.wt)
#column names of the phenotype data
varLabels(dna.e2f.wt)
#expression matrix
exprs(dna.e2f.wt)
#AnnotatedDataFrame
phenoData(dna.e2f.wt)

##subsetting
dna.e2f.wt[1:200,]

##phenoData
dna.e2f.wt <- pheno(dna.e2f.wt,no.groups=2,replicates=c(3,3),groupnames=c("E2F","WT"),DNAreplicates=1)
phenoData(dna.e2f.wt)
pData(dna.e2f.wt)

##useful functions
#vector with pm indices
pmindices<-pmindex(dna.e2f.wt);
#vector with pm sequence
pmseqs<-pmSequence(dna.e2f.wt);
#vector with probe feature names
probens<-probeNames(dna.e2f.wt);
#vector with chromosome info for all PM probes
pmchrs<-pmChr(dna.e2f.wt)

##example wavelet mixed model
pmAllChr<-dna.e2f.wt[pmindices];
#fit - analysis TDDE, no DNA reference
pmAllChrModel<-wmm.fit(pmAllChr,analysis="TDDE",verbose=TRUE,n.levels=12,DNAref=FALSE)
#inference
pmAllChrInference<-wmm.inference(pmAllChrModel,verbose=TRUE);

