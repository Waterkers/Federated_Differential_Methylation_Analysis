if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)

setwd("/home/rstudio") # set the working directory
Betas <- read.csv("/home/rstudio/QC_GSE66351/Normalised_Betas.csv", row.names = 1, header = TRUE)
## change to the relevant input file location
design <- read.csv("/home/rstudio/design_mat/Small_EWAS_design_local.csv", row.names = 1)
identifier = "GSE66351"
## create a folder to save EWAS output
if(!dir.exists(paste0("EWAS_", identifier))){
  dir.create(paste0("EWAS_", identifier))
}
if(!dir.exists(file.path(paste0("EWAS_", identifier), "Plots"))){
  dir.create(file.path(paste0("EWAS_", identifier), "Plots"))
}

model <- lmFit(Beats_csv, design)
contrastMatrics <- makeContrasts(AD - CTRL,levels = colnames(model$coefficients))
contrast_fit <- contrasts.fit(model, contrastMatrics)
result <- eBayes(contrast_fit)
table_res <- topTable(result, adjust="BH",resort.by="P",p.value=1,confint=TRUE,number=dim(Betas)[1])
#save the results to a csv as as
write.csv(table_res, file.path(paste0("EWAS_", identifier),"Small_Results_dataset.csv"))
