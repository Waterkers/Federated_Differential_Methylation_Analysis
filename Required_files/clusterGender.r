#' Predict the sex of samples from DNA methylation data. This function derives pricipal components from the
#' DNA methylation data and uses known sex to identify which PCs separate the sample into two groups that represent sex.
#' It then uses these clusters to predict sex for the whole sample. Therefore it can be used to confirm that sex is correctly #' reported$
#'
#' @param betas A matrix of DNA methylation values.
#' @param sex A vector of sample sexes.
#' @param npcs The number of Principal Components to test.
#' @param thres A numeric threshold to select PCs that correlate with provide sex vector. 
#' @param makePlot A logical to indicate whether a scatterplot of the top two PCs that correlate with sex should be produced.
#' @return A vector of predicted sexes
#' @examples
#' @export
#' @importFrom stats complete.cases cor prcomp kmeans
#' @importFrom graphics plot legend



clusterGender<-function(betas, sex, npcs = 20, thres = 0.5, makePlot = TRUE){
	if(npcs > ncol(betas)){
		npcs<-ncol(betas)
		print(paste("As only", ncol(betas), "samples, can look look at that number of PCs", sep = " "))
	}
	betas.com<-betas[complete.cases(betas),]
	pca<-prcomp(betas.com)

	pca.cor<-rep(NA, npcs)
	for(i in 1:npcs){
		pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)), use = "complete")
	}
	top<-order(abs(pca.cor), decreasing = TRUE)[1]
	second<-order(abs(pca.cor), decreasing = TRUE)[2]
	print(paste("Top correlated principal components with sex:", top, ",", second))
        if(makePlot){
		plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("magenta", "blue")[as.factor(sex)], xlab = paste("PC", top), ylab = paste("PC", second))
		legend("topright", levels(as.factor(sex)), pch = 16, col = c("magenta", "blue"))
	}
	
	predSex_num<-kmeans(pca$rotation[,which(abs(pca.cor) > thres)],2)$cluster
	
	predSex<-rep(NA, length(sex))
	options.sex<-levels(as.factor(sex))
		
	maxSample<-which.max(pca$rotation[,top])
	samplePredSex<-predSex_num[maxSample]
	sampleRepSex<-sex[maxSample]
	sampleOptionSex<-which(options.sex==sampleRepSex)
	
	predSex[which(predSex_num==samplePredSex)]<-options.sex[sampleOptionSex]
	predSex[which(predSex_num!=samplePredSex)]<-options.sex[-sampleOptionSex]
	
	return(predSex)
}

