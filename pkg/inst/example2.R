## example - read.as.FeatureSet ##
##################################

#data.frame PMchr1_raw
#create pData
exprs<-as.matrix(PMchr1_raw[,1:7])
col1<-c(rep("Case",3),rep("Control",3),rep("DNAref",1))
col2<-c(rep("E2F",3),rep("WT",3),rep("DNAref",1))
pData<-data.frame(cbind(col1,col2))
rownames(pData)<-colnames(exprs)
colnames(pData)<-c("Type","Group")
metadata<-data.frame(labelDescription=c("Case/Control Status","Group"),row.names=c("Type","Group"))

tilingExample<-read.as.FeatureSet(exprs=exprs,pkgname="pd.at35b.mr.v04.2.tigrv5",pData=pData,metaData=metadata);


#tilingExample is object of class TilingFeatureSet
class(tilingExample)
#[1] "TilingFeatureSet"
#attr(,"package")
#[1] "oligoClasses"

##Accessing Data Elements
#names of the features
featureNames(tilingExample)
#unique sample identifiers
sampleNames(tilingExample)
#column names of the phenotype data
varLabels(tilingExample)
#expression matrix
exprs(tilingExample)
#AnnotatedDataFrame
phenoData(tilingExample)

##subsetting
tilingExample[1:200,]


##example wavelet mixed model
#fit - analysis TDDE, DNA reference available and used in analysis
pmDNAModel<-wmm.fit(tilingExample,analysis="TDDE",verbose=TRUE,n.levels=12,DNAref=TRUE)
#inference
pmDNAInference<-wmm.inference(pmDNAModel,verbose=TRUE);

