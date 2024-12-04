#' Function to filter a matrix and exclude probes where a common SNP may affect the hybridisation. 
#'
#' @param betas A matrix to filter out probes with snps in the hybridisation sequence. Rownames must be cp identifiers. 
#' @param population Possible options are AFR, AMR, EUR, SAS, EAS which represent African, admixed American, European, South Asian, and East Asian populations. If you are unsure what population or it is a mixture of these then set this parameter to "Unknown"
#' @param maf Minor allelle threshold to apply to SNPs to filter.
#' @return The betas matrix with the relevant rows removed.
#' @examples
#' @export

filterSNPprobesEdit<-function(betas, population = "EUR", maf = 0.05){
	snpProbes <- read.table('/home/silke/Documents/Fed_EWAS/Federated_Differential_Methylation_Analysis/Required_files/mmc1.txt', header = TRUE) #url("https://pmc.ncbi.nlm.nih.gov/articles/instance/4909830/bin/mmc1.txt")
	if(sum(startsWith(rownames(betas), "cg"))/nrow(betas) == 0){
		stop("Rownames are not cg identifiers")
	}
	if(!population %in% c("Unknown", "AFR", "AMR", "EUR", "SAS", "EAS")){
		stop("Specified population not available list")
	} else {
		if(population == "Unknown"){
			snpProbes<-snpProbes[which(snpProbes$AF >= maf & snpProbes$AF <= (1-maf)),]
		} else{
			snpProbes<-snpProbes[which(snpProbes[,paste(population, "AF", sep = "_")] >= maf & snpProbes[,paste(population, "AF", sep = "_")] <= (1-maf)),]
		}
		betas<-betas[!rownames(betas) %in% unique(snpProbes$IlmnID),]
		return(betas)
	}
}

