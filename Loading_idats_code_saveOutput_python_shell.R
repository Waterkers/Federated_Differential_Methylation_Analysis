### get the commands from the subprocess call ###
input <- commandArgs(trailingOnly = TRUE)

## save input commands as local objects
idat <- input[1]
pheno_information <- input[2]
working_dir <- input[3]
identifier <- input[4]

##### Start with installing the required Bioconductor and then packages ########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

need <- c("wateRmelon", "methylumi", "ChAMP")
if (!require(need, quietly = TRUE))
  BiocManager::install(need)
library(methylumi)
library(wateRmelon)


##### source the Exeter functions needed for the pipeline
lapply(list.files("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\Required_files",pattern = "\\.r$",full.names = T),function(x){source(x)})

## Set the working directory
setwd(working_dir) # set working directory

# create the output directory
QC_output <- paste0("QC_", identifier)
if(!dir.exists(QC_output)){
	dir.create(QC_output)
}
if(!dir.exists(file.path(QC_output, "Plots"))){
  dir.create(file.path(QC_output, "Plots"))
} 
print("Output directories created")

#Creating a function around the code to load the data
load_and_save <- function(idat_path, pheno_info) {
# loading in actual data - GSE66351
pheno1 <- read.table(pheno_info)
pheno1 <- t(pheno1)#transpose the imported tabel to the sample characteristics/ids etc are columns and the samples are rows
pheno1 <- as.data.frame(pheno1)
colnames(pheno1)<- pheno1[1,]
pheno1 <- pheno1[2:191,]

print("Phenotype information imported")

Sample_ID <- getBarcodes(idat_path)
pheno1 <- cbind(pheno1, Sample_ID)

pheno1_half <- as.data.frame(pheno1[1:20,])
barcodes_GSE66351_half <- pheno1_half$Sample_ID # my personal laptop cannot deal with all 190 samples so I'm trying it with the first 20 instead

data <- wateRmelon::readEPIC(barcodes = barcodes_GSE66351_half, pdat = pheno1_half, idatPath = idat_path)

save(data, pheno1, file = file.path(QC_output, "methylumiSet_PhenoDatafrma.RData"))
print("Finished reading and saving .idat files")

#save the information stored in the data2 object (the MethylumiSet object thing) into seperate dataframes
setwd(working_dir)

# betas
raw_betas <- betas(data)
write.csv(betas, file.path(QC_output, "Raw_betas.csv"))
# methylated values
raw_methylated <- methylumi::methylated(data)
write.csv(raw_methylated, file.path(QC_output, "Raw_methylated_intensities.csv"))
# unmethylated values
raw_unmethylated <- methylumi::unmethylated(data)
write.csv(raw_unmethylated, file.path(QC_output, "Raw_unmethylated_intensities.csv"))


}

load_and_save(idat, pheno_information)


