### get the commands from the subprocess call ###
input <- commandArgs(trailingOnly = TRUE)


##### Start with installing the required Bioconductor and then packages ########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

need <- c("wateRmelon", "methylumi", "ChAMP")
if (!require(need, quietly = TRUE))
  BiocManager::install(need)
library(methylumi)
library(wateRmelon)

## save input commands as local objects
idat_path <- input[1]
pheno_information <- input[2]
working_dir <- input[3]


##### source the Exeter functions needed for the pipeline
lapply(list.files("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\ExeterEWASPipeline-master\\R",pattern = "\\.r$",full.names = T),function(x){source(x)})
## Set the working directory
setwd(working_dir) # set working directory

# create the output directory
if(!dir.exists("QC_GSE66351_PythonShell")){
  dir.create("QC_GSE66351_PythonShell")
}
if(!dir.exists("QC_GSE66351_PythonShell\\Plots")){
  dir.create("QC_GSE66351_PythonShell\\Plots")
}

print("Output directories created")

# loading in actual data - GSE66351
pheno1 <- read.table(pheno_information)
pheno1 <- t(pheno1)#transpose the imported tabel to the sample characteristics/ids etc are columns and the samples are rows
pheno1 <- as.data.frame(pheno1)
colnames(pheno1)<- pheno1[1,]
pheno1 <- pheno1[2:191,]

write.csv(pheno1, "QC_GSE66351_PythonShell\\pheno_check.csv")
print("Phenotype information imported")

Sample_ID <- getBarcodes(idat_path)
write.csv(Sample_ID, "QC_GSE66351_PythonShell\\barcodes_check1.csv")
pheno1 <- cbind(pheno1, Sample_ID)
pheno1_half <- as.data.frame(pheno1[1:20,])

write.csv(pheno1_half, "QC_GSE66351_PythonShell\\input_check.csv")

barcodes_GSE66351_half <- pheno1_half$Sample_ID # my personal laptop cannot deal with all 190 samples so I'm trying it with the first 20 instead
write.csv(barcodes_GSE66351_half, "QC_GSE66351_PythonShell\\barcodes_check2.csv")
data <- wateRmelon::readEPIC(barcodes = barcodes_GSE66351_half, pdat = pheno1_half, idatPath = idat_path)

write.csv(pheno1_half,"QC_GSE66351_PythonShell\\input_check_postIDAT.csv")
print("idats")

#save the information stored in the data2 object (the MethylumiSet object thing) into seperate dataframes
setwd(working_dir)
# betas
raw_betas <- betas(data)
write.csv(betas, "QC_GSE66351_PythonShell\\Raw_betas.csv")
# methylated values
raw_methylated <- methylumi::methylated(data)
write.csv(raw_methylated, "QC_GSE66351_PythonShell\\Raw_methylated_intensities.csv")
# unmethylated values
raw_unmethylated <- methylumi::unmethylated(data)
write.csv(raw_unmethylated, "QC_GSE66351_PythonShell\\Raw_unmethylated_intensities.csv")

print("Done")




