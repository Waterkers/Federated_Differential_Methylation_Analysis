if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", update = FALSE)
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressPackageStartupMessages(library(limma))

input <- commandArgs(trailingOnly = TRUE)
workingDir <- input[1]
data_dir <-input[2]
design_file <- input[3]
identifier <- input[4]
split <- input[5]
cohort_effect <- input[6]

setwd(workingDir) # set the working directory
#data_dir <- "E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE66351_splits"
#TODO wrap in if statement if split is true otherwise load the normalised beta from data_dir
if(split){
Betas1 <- read.csv(file.path(data_dir, "Split_1_betas.csv"), row.names = 1, header = TRUE)
Betas2 <- read.csv(file.path(data_dir, "Split_2_betas.csv"), row.names = 1, header = TRUE)
Betas3 <- read.csv(file.path(data_dir, "Split_3_betas.csv"), row.names = 1, header = TRUE)
Betas <- cbind(Betas1, Betas2)
Betas <- cbind(Betas, Betas3)
rm(Betas1, Betas2, Betas3)
gc()
  } else {Betas <- read.csv(file.path(data_dir, "Normalised_Betas.csv"), row.names = 1, header = TRUE)}
## change to the relevant input file location
design <- read.csv(design_file, row.names = 1)
design <- design[, !names(design) %in% c("Brain_region")]
#identifier <- "GSE66351"
Betas <- Betas[row.names(design)]
## create a folder to save EWAS output
#TODO wrap in if statement if cohort_effect is true create folder with central otherwise not
# if(cohort_effect) {
#   if (!dir.exists(paste0("EWAS_R_Central_", identifier))) {
#     dir.create(paste0("EWAS_R_Central_", identifier))
#   }
#   if (!dir.exists(file.path(paste0("EWAS_R_Central_", identifier), "Plots"))) {
#     dir.create(file.path(paste0("EWAS_R_Central_", identifier), "Plots"))
#   }
# } else {if (!dir.exists(paste0("EWAS_R_", identifier))) {
#     dir.create(paste0("EWAS_R_", identifier))
#   }
#   if (!dir.exists(file.path(paste0("EWAS_R_", identifier), "Plots"))) {
#     dir.create(file.path(paste0("EWAS_R_", identifier), "Plots"))
#   }
# }
model <- lmFit(Betas, design)
if(cohort_effect) {
  if (split) {
    write.csv(model, file.path(workingDir, "regression_model_central_from_split.csv"))
  } else {write.csv(model, file.path(workingDir, "regression_model_central.csv"))}

}  else {
  if (split) {
    write.csv(model, file.path(workingDir, "regression_model_from_split.csv"))
  } else {write.csv(model, file.path(workingDir, "regression_model.csv"))}
   }
contrastMatrics <- makeContrasts(AD - CTRL,levels = colnames(model$coefficients))
if(cohort_effect) {
  if (split) {
    write.csv(contrastMatrics, file.path(workingDir, "contrast_matrix_central_from_split.csv"))
  } else {write.csv(contrastMatrics, file.path(workingDir, "contrast_matrix_central.csv"))}

} else {
  if (split) {
    write.csv(contrastMatrics, file.path(workingDir, "contrast_matrix_dataset_from_split.csv"))
  } else { write.csv(contrastMatrics, file.path(workingDir, "contrast_matrix_dataset.csv")) } }
contrast_fit <- contrasts.fit(model, contrastMatrics)
if(cohort_effect) {
  if (split) {
    write.csv(contrast_fit, file.path(workingDir, "fitted_contrast_matrix_central_from_split.csv"))
  } else {write.csv(contrast_fit, file.path(workingDir, "fitted_contrast_matrix_central.csv"))}
} else {
  if (split) {
    write.csv(contrast_fit, file.path(workingDir, "fitted_contrast_matrix_dataset_from_split.csv"))
  } else { write.csv(contrast_fit, file.path(workingDir, "fitted_contrast_matrix_dataset.csv")) } }
result <- eBayes(contrast_fit)
if(cohort_effect) {
  if (split) {
  } else {write.csv(result, file.path(workingDir, "eBayes_output_raw_central_from_split.csv"))
}
  } else { if (split) {
  write.csv(result, file.path(workingDir, "eBayes_output_raw_dataset_from_split.csv"))
} else { write.csv(result, file.path(workingDir, "eBayes_output_raw_dataset.csv")) } }
table_res <- topTable(result, adjust="BH",resort.by="P",p.value=1,confint=TRUE,number=dim(Betas)[1])
#save the results to a csv as as
if(cohort_effect) {
  if (split) {
    write.csv(table_res, file.path(workingDir, "Small_Results_dataset_central_from_split.csv"))
  } else {write.csv(table_res, file.path(workingDir, "Small_Results_dataset_central.csv"))}
} else {if (split) {
  write.csv(table_res, file.path(workingDir, "Small_Results_dataset_from_split.csv"))
} else { write.csv(table_res, file.path(workingDir, "Small_Results_dataset.csv")) } }
