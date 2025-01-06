input <- commandArgs(trailingOnly = TRUE)

## save input commands as local objects
idat <- input[1]
pheno_info <- input[2]
working_dir <- input[3]
manifest_path <- input[4]
identifier <- input[5]

##### Start with installing the required packages ########
need <- c("wateRmelon", "methylumi", "ChAMP")
if (!require(need, quietly = TRUE))
  BiocManager::install(need, updata = FALSE)
library(wateRmelon, methylumi)
library(ChAMP) # for some reason library only loads this package when it is called on its own

if (!require("tidyverse", quietly = TRUE))
	install.packages("tidyverse")
library(tidyverse)
##### source the Exeter functions needed for the pipeline - Change to the local filepath that contains these functions
lapply(list.files("/cosybio/project/vanElferen/FedEWAS/Federated_Differential_Methylation_Analysis/Required_files",pattern = "\\.r$",full.names = T),function(x){source(x)})
## Set the working directory
setwd(working_dir) # set working directory
## create a folder for the QC output
QC_output <- paste0("QC_", identifier)
if(!dir.exists(QC_output)){
	dir.create(QC_output)
}
if(!dir.exists(file.path(QC_output, "Plots"))){
  dir.create(file.path(QC_output, "Plots"))
} 

preprocess <- function(idat, pheno_info, intens_threshold) {
# loading in actual data - GSE66351
pheno1 <- read.table(pheno_info)
pheno1 <- t(pheno1)#transpose the imported tabel to the sample characteristics/ids etc are columns and the samples are rows
pheno1 <- as.data.frame(pheno1)
colnames(pheno1)<- pheno1[1,]
pheno1 <- pheno1[2:191,]

print("Phenotype information imported")

Sample_ID <- getBarcodes(idat)

pheno1 <- cbind(pheno1, Sample_ID)
  write.csv(pheno1_half, file.path(QC_output, "Pheno_Info.csv"))
data <- wateRmelon::readEPIC(barcodes = Sample_ID, pdat = pheno1, idatPath = idat)

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

msetEPIC <- data
pheno <- pheno1

########## Start with the QC pipeline from Exeter ##################
### checking methylated and unmethylated intensities #############
### extract sample intensities 
library(methylumi)
m_intensities<-methylumi::methylated(msetEPIC) 
u_intensities<-methylumi::unmethylated(msetEPIC)
## this gives a matrix where each row is a probe and each column a sample

## summarize the intensities of each sample with a single value, the median 
M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

## create a table to store output of QC pipeline
QCmetrics<-cbind(pheno, M.median, U.median) 

# check median deviations +/- 3* median absolute deviation
M_lower_bound <- median(M.median) - 3 * mad(M.median, constant = 1)
M_lower_bound

M_upper_bound <- median(M.median) + 3 * mad(M.median, constant = 1)
M_upper_bound

U_lower_bound <- median(U.median) - 3 * mad(U.median, constant = 1)
U_lower_bound

U_upper_bound <- median(U.median) + 3 * mad(U.median, constant = 1)
U_upper_bound

## 1500 for buccal/saliva. change this to adjust the threshold at which you filter; 2000 for blood
intens.Thres<-intens_threshold

# make PDF histogram of sample intensities
pdf(file.path(QC_output, "Plots", "Sample_Intensity_histogram.pdf"))
par(mfrow = c(1,2))
hist(M.median, xlab = "Median M intensity")
abline(v=M_lower_bound,col=alpha("blue",0.3),lty=2)
abline(v=M_upper_bound,col=alpha("blue",0.3),lty=2)
hist(U.median, xlab = "Median U intensity")
abline(v=U_lower_bound,col=alpha("blue",0.3),lty=2)
abline(v=U_upper_bound,col=alpha("blue",0.3),lty=2)
par(mfrow = c(1,1))
plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity")
abline(v = intens.Thres, col = "red")
abline(h = intens.Thres, col = "red")
abline(v=M_lower_bound,col=alpha("blue",0.3),lty=2)
abline(v=M_upper_bound,col=alpha("blue",0.3),lty=2)
abline(h=U_lower_bound,col=alpha("blue",0.3),lty=2)
abline(h=U_upper_bound,col=alpha("blue",0.3),lty=2)
dev.off()

# State outliers if strict
strictoutliers = names(which(M.median<M_lower_bound | U.median<U_lower_bound | M.median>M_upper_bound | U.median>U_upper_bound))
strictoutliers

# store outliers in QCmetrics
QCmetrics$strict_outliers = QCmetrics$Sample_ID%in%strictoutliers #change column to match the sample_id column name

## calculate a summary statistic for each chip
chip.M.median<-aggregate(M.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(m_intensities), by = 2)]), FUN = median)
chip.U.median<-aggregate(U.median, by = list(unlist(strsplit(colnames(m_intensities), "_"))[seq(from = 1, to = 2*ncol(u_intensities), by = 2)]), FUN = median)



if ("plate" %in% colnames(pheno)) {
	## plot each plate as a boxplot - this does not work for the sample data since there is only a small set
# of samples provided
  pdf(file.path(QC_output, "Plots", "Sample_Intensity_ByPlate_boxplot.pdf"))
  par(mfrow = c(1,2))
  par(mar = c(8, 4, 1, 1))
  nCol<-length(unique(pheno$plate))## assumes there is a column called Plate in your phenotype file
  boxplot(M.median ~ pheno$plate, ylab = "Median M intensity", xlab = "Plate", las = 2, col = rainbow(nCol)) 
  boxplot(U.median ~ pheno$plate, ylab = "Median U intensity", xlab = "Plate", las = 2, col = rainbow(nCol))
  dev.off()
 
## alternatively colour points in original scatterplot by Plate
  pdf(file.path(QC_output, "Plots", "Sample_Intensity_ByPlate_histogram.pdf"))
  nCol<-length(unique(pheno$plate))## assumes there is a column called Plate in your phenotype file
  plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col = rainbow(nCol)[factor(pheno$plate)])
  abline(v = intens.Thres, col = "red")
  abline(h = intens.Thres, col = "red") 
  legend("topright", levels(factor(pheno$plate)), col = rainbow(nCol), pch = 16)
  dev.off() 
	
	 }


##### Checking the Bisulfite conversion ################
bs<-wateRmelon::bscon(msetEPIC) # - doesn't work because of this error "Error in bsI.green[1:2, ] : subscript out of bounds"

pdf(file.path(QC_output, "Plots", "Bisulphite_Conversion.pdf"))
hist(bs, xlab = "Median % BS conversion", main = "")
abline(v = 80, col = "red")
dev.off()
QCmetrics<-cbind(QCmetrics, bs)

##### Checking for biological sex vs. measured sex ############
betas <- methylumi::betas(msetEPIC)
pheno<-pheno[match(colnames(betas), pheno$Sample_ID),]

pdf(file.path(QC_output, "Plots", "Gender_Cluster.pdf"))
predSex1<-findGenderPC(betas, pheno$Sex, npcs = 20) # the default setting of npcs = 20 gave an error so
# I reduced the number of princicple components for the example data set -> remember to put it back at 20
# when using real data
predSex2<-clusterGender(betas, pheno$Sex)
dev.off()

# Confirm findings
PCA = prcomp(betas[complete.cases(betas),])
plot(PCA$rotation[,1],PCA$rotation[,3],col=ifelse(pheno$Sex=="Sex: M","blue","magenta"), main= "PCA of betas, colored by phenotype trait sex",pch=19)
# Check predsex
plot(PCA$rotation[,1],PCA$rotation[,3],col=ifelse(predSex1=="Sex: M","blue","magenta"), main= "PCA of betas, colored by predicted trait predSex1",pch=19)

plot(PCA$rotation[,1],PCA$rotation[,3],col=ifelse(predSex2=="Sex: M","blue","magenta"), main= "PCA of betas, colored by predicted trait predSex2",pch=19)

# the next step would be to put a note in the QCmetrics matrix for the samples where biological/phenotypical 
# sex does not match with the sex determined based on the methylation data - in the example data this is not 
# the case for any samples
# Manual identification - the thresholds are based on the PC plot and visually estimated looking at the clustering
# 2 samples which are classified as male but epigenetically are not. These are located PC3<0 & PC3> -0.1 & PC4 > 0 & PC4 < 0.05
# PC3 = PCA$rotation[,4]
# PC4 = PCA$rotation[,2]
# samples_sex_wrong_shouldbefemale = names(which((PC3 < 0 & PC3 > -0.1) & (PC4 > 0 & PC4 < 0.05) & pheno$sexem1=="MALE"))
# samples_sex_wrong_shouldbefemale
# 2 samples which are classified as female but epigenetically are not. These are located PC3 > 0 & PC3 < 0.1 & PC4 > -0.15 & PC4 < -0.05
# samples_sex_wrong_shouldbemale= names(which((PC3 > 0 & PC3 < 0.1 ) & (PC4 > -0.15 & PC4 < -0.05) & pheno$sexem1=="FEMALE"))
# samples_sex_wrong_shouldbemale

# Saving the correct predicted sex results to the QCmetrics matrix
QCmetrics<-cbind(QCmetrics, predSex1)

###### Check genetically identical samples correlate across SNP probes ######
#- Does not work for the example data
betas <- methylumi::betas(msetEPIC) 
pheno<-pheno[match(colnames(betas), pheno$SampleLabel),]

betas.rs<-betas[grep("rs", rownames(betas)),]

# check for complete cases as it will cause error
betas.rs = betas.rs[complete.cases((betas.rs)),]

snpCor<-cor(betas.rs)
#names(snpCor)<-pheno$sentrix_full ## the deafult here is the sentrix_full of the sample which you may wish to change to a more user friendly identifier

# remove correlation to itself
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}

# calculate the maximum correlation for each sample with all other samples (except for itself)
corMax<-apply(snpCor, 1, max, na.rm = TRUE)

pdf(file.path(QC_output, "Plots", "SNP_Correlations.pdf"))
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()

# Check which samples match to their best match 
pdf(file.path(QC_output, "Plots", "SNP_Samples.pdf"), width = 15, height = 8)
par(mfrow = c(2,4))
for(i in 1:ncol(betas.rs)){
  val = betas.rs[,i]
  cors = cor(val,betas.rs)
  cors[i] = NA
  o = which.max(cors)
  xname = paste0(colnames(betas.rs)[i]," (ID: ",pheno$Sample_title[which(pheno$Sample_ID==colnames(betas.rs)[i])],")")
  yname = paste0( colnames(cors)[o]," (ID: ",pheno$Sample_title[which(pheno$Sample_ID==colnames(cors)[o])],")")
  
  id_match = pheno$Sample_title[which(pheno$Sample_ID==colnames(betas.rs)[i])] == pheno$Sample_title[which(pheno$Sample_ID==colnames(cors)[o])]
  
  id_match_text = ifelse(id_match,"IDs match","IDs dont match")
  #if(!id_match){warning(paste0("mismatch with sample ",xname, " Index: ",i))}
  
  plot(betas.rs[,i], betas.rs[,o],xlab=xname,ylab=yname,main=paste0(id_match_text,"\n Cor: ",round(cors[o],3)))
  
}
dev.off()

# Storing the correlations, source, target
duplicateSamples = data.frame(source = rownames(snpCor), target = names(apply(snpCor, 1, max, na.rm = TRUE)), maxcorrelation = as.numeric(apply(snpCor, 1, max, na.rm = TRUE)))


QCmetrics<-cbind(QCmetrics, duplicateSamples)

##### P-filter #############
# Remove objects in memory as a lot of RAM is needed when running pfilter
removelist = c("betas","chip.M.median","chip.U.median","cors","duplicateSamples","m_intensities","PCA","u_intensities","snpCor","betas.rs","bs","corMax","i","id_match","id_match_text","M.median","o","predSex1","predSex2","U.median","val","xname","yname")
rm(list=removelist)
gc()

# filter on bad samples and sites (CpGs)
pdf(file.path(QC_output, "Plots", "Betas_raw_boxplot.pdf"))
boxplot(methylumi::betas(msetEPIC),main="Betas raw betas")
dev.off()


# Filter on quality
msetEPIC.pf <- wateRmelon::pfilter(msetEPIC, perc = 5) # samples to 5% due to buccal
# reload methylumi package when running into issues
retained_probes <- rownames(betas(msetEPIC)) %in% rownames(betas(msetEPIC.pf))
# store info
pFilterPass<-colnames(betas(msetEPIC)) %in% colnames(betas(msetEPIC.pf))
QCmetrics<-cbind(QCmetrics, pFilterPass)

# make room in ram
rm(msetEPIC)
gc()

##### Removal of cross-hybridisation probes ###############
# load the cpgs to be removed - this step should probably be preformed locally
crosshyb <- read.table("/cosybio/project/vanElferen/FedEWAS/Federated_Differential_Methylation_Analysis/Required_files/mmc2.txt")
snpProbes <- read.table("/cosybio/project/vanElferen/FedEWAS/Federated_Differential_Methylation_Analysis/Required_files/mmc1.txt", header = TRUE)
msetEPIC.pf <- msetEPIC.pf[!(rownames(msetEPIC.pf@assayData$betas) %in% crosshyb[,1]), ]
kept_probes <- filterSNPprobesEdit(msetEPIC.pf, population = "EUR", maf = 0.05)
msetEPIC.pf <- msetEPIC.pf[rownames(msetEPIC.pf@assayData$betas) %in% rownames(kept_probes), ]
msetEPIC.pf <- msetEPIC.pf[-grep("rs", rownames(msetEPIC.pf@assayData$betas)), ]

#save the information stored in the data2 object (the MethylumiSet object thing) into seperate dataframes
setwd(working_dir)

# save the filtered MethylumiSet object for r-based normalisation
save(msetEPIC.pf, QCmetrics, file = file.path(QC_output,"Filtered_MethyLumiSet.RData"))

# betas
raw_betas <- betas(msetEPIC.pf)
write.csv(raw_betas, file.path(QC_output, "Filtered_Betas.csv"))
# methylated values
raw_methylated <- methylumi::methylated(msetEPIC.pf)
write.csv(raw_methylated, file.path(QC_output, "Filtered_Methylated.csv"))
# unmethylated values
raw_unmethylated <- methylumi::unmethylated(msetEPIC.pf)
write.csv(raw_unmethylated, file.path(QC_output, "Filtered_Unmethylated.csv"))

write.csv(QCmetrics, file.path(QC_output, "pre_norm_pheno_information.csv"))


}

preprocess(idat, pheno_info, 2000)


