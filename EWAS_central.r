if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)

setwd("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\Data_Full_Datasets\\GSE66351\\EWAS") # set the working directory
data_dir <- "E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE66351_splits"
Betas1 <- read.csv(file.path(data_dir, "Split_1_betas.csv"), row.names = 1, header = TRUE)
Betas2 <- read.csv(file.path(data_dir, "Split_2_betas.csv"), row.names = 1, header = TRUE)
Betas3 <- read.csv(file.path(data_dir, "Split_3_betas.csv"), row.names = 1, header = TRUE)
Betas <- cbind(Betas1, Betas2)
Betas <- cbind(Betas, Betas3)
rm(Betas1, Betas2, Betas3)
gc()
## change to the relevant input file location
design <- read.csv("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE66351_splits\\central_design_matrix.csv", row.names = 1)
design <- design[, !names(design) %in% c("Brain_region")]
identifier = "GSE66351"
Betas <- Betas[row.names(design)]
## create a folder to save EWAS output
if(!dir.exists(paste0("EWAS_", identifier))){
  dir.create(paste0("EWAS_", identifier))
}
if(!dir.exists(file.path(paste0("EWAS_", identifier), "Plots"))){
  dir.create(file.path(paste0("EWAS_", identifier), "Plots"))
}

model <- lmFit(Betas, design)
contrastMatrics <- makeContrasts(AD - CTRL,levels = colnames(model$coefficients))
contrast_fit <- contrasts.fit(model, contrastMatrics)
result <- eBayes(contrast_fit)
table_res <- topTable(result, adjust="BH",resort.by="P",p.value=1,confint=TRUE,number=dim(Betas)[1])
#save the results to a csv as as
write.csv(table_res, file.path(paste0("EWAS_", identifier),"Small_local_Results_dataset_new.csv"))
