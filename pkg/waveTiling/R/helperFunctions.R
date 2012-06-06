
## functions for data input

cel2TilingFeatureSet <- function(dataPath,annotationPackage)
{
	expFiles <- paste(dataPath,list.celfiles(dataPath),sep="/")
	tfs <- read.celfiles(expFiles,pkgname=annotationPackage)
	return(tfs)
}


extractChrFromFasta <- function(char,fasPath)
{
	splt <- unlist(strsplit(char,fasPath))
	hlp <- splt[seq(2,length(splt),2)]
	splt2 <- unlist(strsplit(hlp,".fas"))
	return(splt2)
}

pm2mm <- function(z)
{
	z <- complement(subseq(z,start=13,end=13))
	return(z)
}

## functions for model fitting

normVec <- function(vec)
{
	return(vec/sqrt(sum(vec^2)))
}

wave.transform <- function(x,n.levels=log(length(x),2),filt="haar") 
{
	N <- length(x)
	J <- n.levels
	if (N/2^J != trunc(N/2^J)) stop("Sample size is not divisible by 2^J")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf
	storage.mode(h) <- "double"
	g <- dict$lpf
	storage.mode(g) <- "double"
	start <- 1
	stop <- N/2
	Y <- rep(0,N)
	for (j in 1:J)
	{
		W <- V <- numeric(N/2^j)
		out <- .C("dwt", as.double(x), as.integer(N/2^(j - 1)),L, h, g, W =  as.double(W), V = as.double(V), PACKAGE = "waveslim")[6:7]
		Y[start:stop] <- out$W
		x <- out$V
		start <- start+N/(2^j)
		stop <- stop+N/2^(j+1)
	}
	stop <- N
	Y[start:stop] <- x
	return(Y)
}

wave.backtransform <- function(D,n.levels,filt="haar")
{
	p <- length(D)
	if ((p%%2^n.levels)!=0) stop("Sample size is not divisible by 2^n.levels")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf
	storage.mode(h) <- "double"
	g <- dict$lpf
	storage.mode(g) <- "double"
	stop <- p
	start <- p-p/2^n.levels+1
	X <- D[start:stop]
	for (j in n.levels:1)
	{
		stop <- start-1
		start <- stop-p/2^j+1
		N <- length(X)
		XX <- numeric(2*p/2^j)
		X <- .C("idwt", as.double(D[start:stop]), as.double(X), as.integer(N),L, h, g, out = as.double(XX), PACKAGE = "waveslim")$out
	}
	return(X)
}

wave.backtransformK <- function(Cum,n.levels,order=1,filt="haar")
{
	p <- length(Cum)
	if ((p%%2^n.levels)!=0) stop("Sample size is not divisible by 2^n.levels")
	dict <- wave.filter(filt)
	L <- dict$length
	storage.mode(L) <- "integer"
	h <- dict$hpf^order
	storage.mode(h) <- "double"
	g <- dict$lpf^order
	storage.mode(g) <- "double"
	stop <- p
	start <- p-p/2^n.levels+1
	X <- Cum[start:stop]
	for (j in n.levels:1)
	{
		stop <- start-1
		start <- stop-p/2^j+1
		N <- length(X)
		XX <- numeric(2*p/2^j)
		X <- .C("idwt", as.double(Cum[start:stop]), as.double(X), as.integer(N),L, h, g, out = as.double(XX), PACKAGE = "waveslim")$out
	}
	return(X)
}

MapMar <- function(D,X,vareps,J=0,eqsmooth=TRUE)
{
	N <- nrow(X)
	q <- ncol(X)
	K <- ncol(D)
	if (J==0)
	{
		ends <- K
	} else
	{
	ends <- cumsum(K/2^c(1:J,J))
	}
	B <- matrix(0,nrow=q,ncol=K)
	varB <- B
	phi <- B
	if (eqsmooth)
	{
		out <- .C("MAPMARGEQSMOOTH",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	} else
	{
		out <- .C("MAPMARG",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	}
	names(out) <- c("beta_MAP","varbeta_MAP","phi")
	dim(out[[1]]) <- c(K,q)
	out[[1]] <- t(out[[1]])
	dim(out[[2]]) <- c(K,q)
	out[[2]] <- t(out[[2]])
	dim(out[[3]]) <- c(K,q)
	out[[3]] <- t(out[[3]])
	return(out)
}

MapMarImp <- function(D,X,vareps,J=0,eqsmooth=TRUE)
{
	N <- nrow(X)
	q <- ncol(X)
	K <- ncol(D)
	if (J==0)
	{
		ends <- K
	} else
	{
	ends <- cumsum(K/2^c(1:J,J))
	}
	B <- matrix(0,nrow=q,ncol=K)
	varB <- B
	phi <- B
	if (eqsmooth)
	{
		out <- .C("MAPMARGIMPEQSMOOTH",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	} else
	{
		out <- .C("MAPMARGIMP",as.double(D),as.integer(K),as.double(vareps),as.double(B),as.double(varB),as.double(phi),as.double(X),as.double(diag(t(X)%*%X)),as.integer(q),as.integer(N),as.integer(ends), PACKAGE = "waveTiling")[4:6]
	}
	names(out) <- c("beta_MAP","varbeta_MAP","phi")
	dim(out[[1]]) <- c(K,q)
	out[[1]] <- t(out[[1]])
	dim(out[[2]]) <- c(K,q)
	out[[2]] <- t(out[[2]])
	dim(out[[3]]) <- c(K,q)
	out[[3]] <- t(out[[3]])
	return(out)
}


WaveMarEstVarJ <- function(Y,X,n.levels,wave.filt,prior=c("normal","improper"),eqsmooth=TRUE,saveall=FALSE,D,var.eps,max.it,tol=1e-6,trace=FALSE)
{
	N <- nrow(D)
	K <- ncol(D)		
	Kj <- K/2^c(1:n.levels,n.levels)
	ends <- cumsum(Kj)
	starts <- c(1,ends[-(n.levels+1)]+1)
	varEps <- rep(mad(D[,starts[1]:ends[1]])^2,n.levels+1)
	if (var.eps=="mad")
	{
		max.it <- -1
		saveall <- TRUE
	}
	if (max.it<0)
	{
		for (i in 2:(n.levels+1))
		{
			varEps[i] <- ifelse(i<(n.levels+1),mad(D[,starts[i]:ends[i]])^2,varEps[i-1])
		}
	}
	crit1 <- sum(D^2)
	if (max.it<1)
	{
		if (prior=="normal")
		{
			WaveFit <- MapMar(D,X,varEps,n.levels,eqsmooth)
		}
		if (prior=="improper")
		{
			WaveFit <- MapMarImp(D,X,varEps,n.levels,eqsmooth)
		}
	} else
	{
		if (trace==TRUE) message("Gauss Seidel algorithm...")
		for (j in 1:max.it)
		{
			crit0 <- crit1
			if (prior=="normal")
			{
				WaveFit <- MapMar(D,X,varEps,n.levels,eqsmooth)
			}
			if (prior=="improper")
			{
				WaveFit <- MapMarImp(D,X,varEps,n.levels,eqsmooth)
			}
			for (i in 1:(n.levels+1))
			{
				VinvD <- D[,starts[i]:ends[i]]-X%*%(t(X)%*%D[,starts[i]:ends[i]]/(diag(t(X)%*%X)+1/WaveFit$phi[,starts[i]:ends[i]]))
				varEps[i] <- sum(D[,starts[i]:ends[i]]*VinvD)/(ends[i]-starts[i]+1)/N
			}
			crit1 <- -1/2*sum(log(t(X)%*%X%*%WaveFit$phi+rep(1,ncol(X))))-1/2*(log(varEps)%*%Kj)-N*K*log(2*pi)/2-N*K/2
			if (trace==TRUE) message("iteration ", j," of ",max.it)
			if ((abs((crit0-crit1)/crit0))<tol)
			{
				if (trace==TRUE) message("Gauss Seidel algorithm converged")
				break
			}
		}
		if ((abs(crit0-crit1)/crit0)>tol) message("Warning: The maximum number of iterations reached without convergence")
	}
	if (!saveall)
	{
		WaveFit$beta_MAP <- matrix()
		WaveFit$varbeta_MAP <- matrix()
		WaveFit$phi <- matrix()
	}
	WaveFit$varEps <- varEps
	WaveFit$n.levels <- n.levels
	WaveFit$K <- K
	WaveFit$N <- N
	WaveFit$Kj <- Kj
	WaveFit$X <- X
	WaveFit$wave.filt <- wave.filt
	return(WaveFit)
}

## plot

makeNewAnnotationTrack <- function(biomartObj,chromosome,strand,minBase,maxBase,feature="exon",dp=NULL)
{
	if (is.null(dp))
	{
		dp <- DisplayPars(exon="lightblue")
	}
	exn <- exons(biomartObj)
	exnId <- which((as.character(seqnames(exn))==as.character(chromosome)&(as.character(strand(exn))==strand)&(start(ranges(exn))>minBase)&(end(ranges(exn))<maxBase)))
	exnSel <- exn[exnId]
	regionObj <- data.frame(start=start(ranges(exnSel)),end=end(ranges(exnSel)),feature=rep(feature,length(exnSel)),group=elementMetadata(exnSel)$exon_id,ID=rep(NA,length(exnSel)))
	return(makeAnnotationTrack(regions=regionObj,dp=dp))
}


makeNewAnnotationTextOverlay <- function(biomartObj,chromosome,strand,minBase,maxBase,region,feature=c("gene","transposable_element_gene"),y=.5,dp=NULL)
{
	if (is.null(dp))
	{
		dp=DisplayPars(cex=1)
	}
	trs <- transcripts(biomartObj)
	trsId <- which((as.character(seqnames(trs))==as.character(chromosome)&(as.character(strand(trs))==strand)&(start(ranges(trs))>minBase)&(end(ranges(trs))<maxBase)))
	trsSel <- trs[trsId]
	annohlp <- data.frame(ID=elementMetadata(trsSel)$tx_name,start=start(ranges(trsSel)),end=end(ranges(trsSel)))
	return(makeTextOverlay(as.character(annohlp$ID),xpos=(annohlp$start+annohlp$end)/2,ypos=y,region=region,dp=dp))
}

makeNewTranscriptRectangleOverlay <- function(sigRegions,locations,start,end,region=NULL,dp=NULL)
{
	if (is.null(dp))
	{
		dp=DisplayPars(color="black")
	}
	detectedRegionsSelect <- matrix(sigRegions[sigRegions[,1]<end&sigRegions[,2]>start,],ncol=2)
	detectedRegionsSelect[detectedRegionsSelect[,1]<start,1] <- start
	detectedRegionsSelect[detectedRegionsSelect[,2]>end,2] <- end
	return(makeRectangleOverlay(start=locations[detectedRegionsSelect[,1]],end =locations[detectedRegionsSelect[,2]],region=region,dp=dp)) 
}

#### makeDesign
makeDesign<-function(design=c("time","circadian","group","factorial"),replics, noGroups, factor.levels=NULL) {
# 	if (method=="twoGroup")
# 	{
# 		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=2)
# 		Xorig[,1] <- 1
# 		Xorig[,2] <- rep(c(1,-1),replics) ## contr.helmert
# 	}
	if (design=="time") #time
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=noGroups)
		orderedFactor <- factor(1:noGroups,ordered=TRUE)
		desPoly <- lm(rnorm(noGroups)~orderedFactor,x=TRUE)$x
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desPoly[,2:noGroups],2,rep,replics)
		
	}
	else if (design=="circadian") # circadian
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=3)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=3)
		Xorig[,1] <- 1
		Xorig[,2] <- rep(sin(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
		Xorig[,3] <- rep(cos(seq(0,2*pi-(pi/noGroups),2*pi/noGroups)),replics)
	}
	else  if (design=="group") # group
	{
#		Xorig <- matrix(0,nrow=dim(pData(object))[1],ncol=noGroups)
		Xorig <- matrix(0,nrow=noGroups*replics,ncol=noGroups)
		desHelmert <- contr.helmert(noGroups)
		Xorig[,1] <- 1
		Xorig[,2:noGroups] <- apply(desHelmert[,1:(noGroups-1)],2,rep,replics)
	}
	else if (design=="factorial")
	{
		if (is.null(factor.levels))
		{
			stop("No proper factor levels given.")
		}
		Xorig <- matrix(0,nrow=sum(replics),ncol=prod(factor.levels))
		Xorig[,1] <- 1
		desTrt1 <- contr.treatment(factor.levels[1])
		desTrt2 <- contr.treatment(factor.levels[2])
		#desHelmert1 <- contr.helmert(factor.levels[1])
		#desHelmert2 <- contr.helmert(factor.levels[2])
		replicsMat <- matrix(replics,nrow=factor.levels[1],ncol=factor.levels[2],byrow=TRUE)
		desFac1 <- apply(desTrt1,2,function(x) unlist(sapply(1:nrow(replicsMat),function(y) rep(x[y],sum(replicsMat[y,])))))
		Xorig[,2:(factor.levels[1])] <- desFac1
		desFac2 <- apply(desTrt2,2,function(x) unlist(sapply(1:nrow(replicsMat),function(y) rep(x,replicsMat[y,]))))
		Xorig[,(factor.levels[1]+1):(factor.levels[1]+factor.levels[2]-1)] <- desFac2
		# include interactions
		colno <- sum(factor.levels)
		for (i in 1:ncol(desFac1))
		{
			for (j in 1:ncol(desFac2))
			{
				Xorig[,colno] <- desFac1[,i]*desFac2[,j]
				colno <- colno + 1
			}
		}
	}
	else {
		stop("No proper design specified!")
	}	
	out <- list()
	X.qr <- qr(Xorig)
	out$Xorthnorm <- qr.Q(X.qr)
	out$Xorig <- Xorig
	return (out);
}

#### makeContrasts
makeContrasts <- function(contrasts, nlevels) {
	if (contrasts=="compare") {
		q <- nlevels*(nlevels-1)/2
		contr <- matrix(0,nrow=q,ncol=nlevels)
		hlp1 <- rep(2:nlevels,1:(nlevels-1))
		hlp2 <- unlist(sapply(1:(nlevels-1),function(x) seq(1:x)))
		for (i in 1:nrow(contr))
		{
			contr[i,hlp1[i]] <- 1
			contr[i,hlp2[i]] <- -1
		}
	}

	return (contr);
}


