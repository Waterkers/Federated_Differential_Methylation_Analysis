############################################################################
#                         EWAS                                             #
###########################################################################
# load the exeter functions
lapply(list.files("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\ExeterEWASPipeline-master\\R",pattern = "\\.r$",full.names = T),function(x){source(x)})
# setwd("") # set the working directory
load("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\R_Pipeline\\QC\\GSE66351_first20_QCandDasen.RData")
## create a folder to save EWAS output
if(!dir.exists("EWAS")){
  dir.create("EWAS")
}
if(!dir.exists("EWAS/Plots")){
  dir.create("EWAS/Plots")
}

#### Get rid of cross-hybridising probes, based on McCartney et al., 2016 ###
crosshyb <- read.table("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\Cross_hybridising_CpGTargetting_Probes_McCartneyetal2016.txt")

betas<-Betas[!(rownames(Betas) %in% crosshyb[,1]), ]
betas<-filterSNPprobes(betas, population = "EUR", maf = 0.05) ## filters common probes based on allele frequency in european populations.
betas<-betas[-grep("rs", rownames(betas)),] ## remove SNP probes

##### Run the EWAS ##########
res<-matrix(data = NA, nrow = nrow(betas), ncol = 3)
colnames(res)<-c("Diagnosis_Beta", "Diagnosis_SE", "Diagnosis_P")
rownames(res)<-rownames(betas)
# removed the smoking score variable from the linear model because not all data was provided to run the 
# smokingScore function and calculate the smoking score for the samples. 
for(i in 1:nrow(betas)){
  model<-lm(betas[i,] ~ Small_Pheno$Diagnosis + Small_Pheno$Age + factor(Small_Pheno$Sex) + Small_Pheno$Cell_Type.CT1 + Small_Pheno$Cell_Type.CT2 + Small_Pheno$Cell_Type.CT3 + factor(Small_Pheno$Sentrix_ID))
  res[i,c(1)]<-coefficients(model)["Small_Pheno$Diagnosis"]
  res[i,2]<-summary(model)$coefficients["Small_Pheno$Diagnosis",2]
  res[i,c(3)]<-summary(model)$coefficients["Small_Pheno$Diagnosis",4]
}

# the model part of the code runs, but it doesn't provide a standard error or p-value for each row
# meaning that the loop doesn't work and there is no results matrix generated