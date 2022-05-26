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
library(limma)

setwd("/home/rstudio") # set the working directory
load("/home/rstudio/QC_GSE66351/Preprocessed_CTD.RData")
load("/home/rstudio/QC_GSE66351/Full_Phenotype_Information.RData")# change to the relevant input file location
design <- read.csv("/home/rstudio/design_mat/Small_EWAS_design_local.csv", row.names = 1)
identifier = "GSE66351"
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

oc_design <- read.csv("design_mat/Small_Occipitalcortex_EWAS_design_local.csv", row.names = 1)
ft_design <- read.csv("design_mat/Small_Frontalcortex_EWAS_design_local.csv", row.names = 1)
tc_design <- read.csv("design_mat/Small_Temporalcortex_EWAS_design_local.csv", row.names = 1)
tissues <- unique(Full_Pheno[tissue_column])
tissue_design <- list(oc_design, ft_design, tc_design)
for (i in 1:dim(tissues)[1]) {
  Betas_t <- Betas[ ,Full_Pheno[tissue_column] == tissues[i,1]]
  design_t <- tissue_design[i]
  
##### Run the EWAS ##########

  # limma version of linear model
  
    model <- lmFit(Betas_t, design_t[[1]])
    contrastMatrics <- makeContrasts(AD - CTRL,levels = colnames(model$coefficients))
    contrast_fit <- contrasts.fit(model, contrastMatrics)
    result <- eBayes(contrast_fit)
    table_res <- topTable(result, adjust="BH",resort.by="P",p.value=1,confint=TRUE,number=dim(Betas)[1])
#save the results to a csv as as
    write.csv(table_res, file.path(paste0("EWAS_", identifier), paste0(tissues[i,1], "_","Results_dataset.csv")))

#add annotation
  data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  data("Locations")
  force(Locations)
  included_annotations <- cbind(c(Locations["chr"], Islands.UCSC, Other[c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "HMM_Island")]))
  included_annotations <- included_annotations[match(rownames(Betas_t), rownames(included_annotations)),]
  res_annotated <- cbind(table_res, included_annotations)


#save the results to a csv as as
  write.csv(res_annotated, file.path(paste0("EWAS_", identifier), paste0(tissues[i,1], "_", "Annotated_Results_dataset.csv")))
}
#save the results in the standard format defined for a BED file - needed for the DMR calling function
# standard format that requires the first three columns to be: chrom, start and end. 
# This information was downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/
# for the method the 4th column is assumed to contain the p-value for a region
  standard_bed_columns <- read.table("/home/rstudio/HAIB.A549.EtOH.Rep.3.bed")
  standard_bed_columns <- standard_bed_columns[ ,c(1,2,3,4)]

  res_bed <- merge(standard_bed_columns, table_res, by.x = 4, by.y = "row.names", all.y = TRUE)

  res_bed <- data.frame(chrom = res_bed[2], start = res_bed[3], stop = res_bed[4], p_value = res_bed$Diagnosis_P, coeffi = res_bed$Diagnosis_Beta, stan_er = res_bed$Diagnosis_SE, Illumina_ID = res_bed[1])
  write.table(res_bed, file.path(paste0("EWAS_", identifier), "results.bed"))

