#' Predict the sex of samples from DNA methylation data. This function derives pricipal components from the 
#' DNA methylation data and uses known sex to identify the top PC associated with sex. It when uses this PC to predict sex in the whole sample
#' The prediction is done by assuming that the PC is centred on 0 and males and females lie either side of this central point.
#'
#'
#' @param betas A matrix of DNA methylation values. 
#' @param sex A vector of sample sexes.
#' @param npcs The number of Principal Components to test.
#' @return A vector of predicted sexes
#' @examples
#' @export
#' @importFrom stats cor
#' @importFrom graphics plot legend


findGenderPC<-function(betas, sex, npcs = 20){

	betas.com<-betas[complete.cases(betas),]
	pca<-prcomp(betas.com)

	pca.cor<-rep(NA, npcs)
	for(i in 1:npcs){
		pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)), use = "complete")
	}
	top<-order(abs(pca.cor), decreasing = TRUE)[1]
	second<-order(abs(pca.cor), decreasing = TRUE)[2]
	print(paste("Top correlated principal components with sex:", top, ",", second))
	plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("magenta", "blue")[as.factor(sex)], xlab = paste("PC", top), ylab = paste("PC", second))
	legend("topright", levels(as.factor(sex)), pch = 16, col = c("magenta", "blue"))
	predSex<-rep(NA, length(sex))
	options.sex<-levels(as.factor(sex))
	
	if(abs(pca.cor[top]) > 0.9){
		print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes")
	} else {
	  print(paste("Top PC has r =", round(abs(pca.cor[top]),2), "with sex so may not be good enough to confirm reported sexes"))
	}
	
	if(sign(pca.cor[top]) == 1){
		predSex[which(pca$rotation[,top] < 0)]<-options.sex[1]
		predSex[which(pca$rotation[,top] > 0)]<-options.sex[2]
	} else {
		predSex[which(pca$rotation[,top] < 0)]<-options.sex[2]
		predSex[which(pca$rotation[,top] > 0)]<-options.sex[1]
	}
	
	return(predSex)
}

