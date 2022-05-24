#### keep in mind before running #######
# change the variable names of the variables in the linear model
# to the desired ones present in the phenotype file that is being used

# change the column to be selected from the model summary object so
# it matches with the system input used in the linear model.
# The general format for a factor variable - diagnosis in this case is:
#		variable_name first_level_factor (Diagnosis CTRL)
# The general format for variables that are converted to numerical variables inline - as defined here is:
#		as.numeric(variable_name)
# The general format for a predefined numerical variable is:
#		variable_name

############################################################################
#                         EWAS                                             #
###########################################################################
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


setwd("/home/rstudio") # set the working directory
load("/home/rstudio/QC_GSE105109/GSE105109_Preprocessed_CTD.RData") # change to the relevant input file location
identifier = "GSE105109"
## create a folder to save EWAS output
if(!dir.exists(paste0("EWAS_", identifier))){
  dir.create(paste0("EWAS_", identifier))
}
if(!dir.exists(file.path(paste0("EWAS_", identifier), "Plots"))){
  dir.create(file.path(paste0("EWAS_", identifier), "Plots"))
}

# split datasets based on source tissue
if (identifier == "GSE66351") {
  tissue_column = "Brain_region"
} else {
  tissue_column = "Source_tissue"
}
tissues <- unique(Betas[tissue_column])
for (i in 1:length(tissues)) {
  Betas <- Betas[Betas[tissue_column] == tissues[i], ]

  
##### Run the EWAS ##########
  res<-matrix(data = NA, nrow = nrow(Betas), ncol = 3)
  colnames(res)<-c("Diagnosis_Beta", "Diagnosis_SE", "Diagnosis_P")
  rownames(res)<-rownames(Betas)

  for(i in 1:nrow(Betas)){
    model<-lm(Betas[i,] ~ Small_Pheno$Diagnosis + as.numeric(Small_Pheno$Age) + factor(Small_Pheno$Sex) + as.numeric(Small_Pheno$Cell_Type.CT1)+ as.numeric(Small_Pheno$Cell_Type.CT2) + as.numeric(Small_Pheno$Cell_Type.CT3)+ factor(Small_Pheno$Sentrix_ID))
    res[i,c(1)]<-coefficients(model)["Small_Pheno$Diagnosis Control"]
    res[i,2]<-summary(model)$coefficients["Small_Pheno$Diagnosis Control",2]
    res[i,c(3)]<-summary(model)$coefficients["Small_Pheno$Diagnosis Control",4]
  }

# adding a column  for with the corrected p-values (Benjanimi-Hochberg)
  corr_pval <- p.adjust(c(res[ ,3]), method = "BH")
  res <- cbind(res, corr_pval)

#save the results to a csv as as
  write.csv(res, file.path(paste0("EWAS_", identifier, "_", tissues[i]), "Results_dataset.csv"))

#add annotation
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  data("Locations")
  force(Locations)
  included_annotations <- cbind(c(Locations["chr"], Islands.UCSC, Other[c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "HMM_Island")]))
  included_annotations <- included_annotations[match(rownames(Betas), rownames(included_annotations)),]
  res_annotated <- cbind(res, included_annotations)


#save the results to a csv as as
  write.csv(res_annotated, file.path(paste0("EWAS_", identifier, "_", tissues[i]), "Annotated_Results_dataset.csv"))

#save the results in the standard format defined for a BED file - needed for the DMR calling function
# standard format that requires the first three columns to be: chrom, start and end. 
# This information was downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/
# for the method the 4th column is assumed to contain the p-value for a region
  standard_bed_columns <- read.table("/home/rstudio/HAIB.A549.EtOH.Rep.3.bed")
  standard_bed_columns <- standard_bed_columns[ ,c(1,2,3,4)]

  res_bed <- merge(standard_bed_columns, res, by.x = 4, by.y = "row.names", all.y = TRUE)

  res_bed <- data.frame(chrom = res_bed[2], start = res_bed[3], stop = res_bed[4], p_value = res_bed$Diagnosis_P, coeffi = res_bed$Diagnosis_Beta, stan_er = res_bed$Diagnosis_SE, Illumina_ID = res_bed[1])
  write.table(res_bed, file.path(paste0("EWAS_", identifier, "_", tissues[i]), "results.bed"))
}
