############################################################################
#                         EWAS                                             #
###########################################################################
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# load the exeter functions
lapply(list.files("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\ExeterEWASPipeline-master\\R",pattern = "\\.r$",full.names = T),function(x){source(x)})
# setwd(" E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data") # set the working directory
load("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\Output\\QC_GSE66351_PythonShell\\Preprocessed_Normalised_MethyLumiSet.RData")
## create a folder to save EWAS output
if(!dir.exists("EWAS_GSE66351")){
  dir.create("EWAS_GSE66351")
}
if(!dir.exists("EWAS_GSE66351/Plots")){
  dir.create("EWAS_GSE66351/Plots")
}

#### Get rid of cross-hybridising probes, based on McCartney et al., 2016 ###
crosshyb <- read.table("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\Cross_hybridising_CpGTargetting_Probes_McCartneyetal2016.txt")
snpProbes <- read.table("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\mmc1.txt", header = TRUE)
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
  model<-lm(betas[i,] ~ factor(Small_Pheno$Diagnosis) + as.numeric(Small_Pheno$Age) + factor(Small_Pheno$Sex) + factor(Small_Pheno$Sentrix_ID))
  res[i,c(1)]<-coefficients(model)["factor(Small_Pheno$Diagnosis) CTRL"]
  res[i,2]<-summary(model)$coefficients["factor(Small_Pheno$Diagnosis) CTRL",2]
  res[i,c(3)]<-summary(model)$coefficients["factor(Small_Pheno$Diagnosis) CTRL",4]
}

# By changing the name of the coefficient to select to its actual name in the table everything now works fine
# However, the fact that I had to manually copy and paste the name from the coefficients names table probably means
# it cannot be easilly automated, at least not in a way that I can think of

#save the results to a csv as as
write.csv(res, "EWAS_GSE66351/Results_dataset_test.csv")

#add annotation
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
data("Locations")
force(Locations)
included_annotations <- cbind(c(Locations["chr"], Islands.UCSC, Other[c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "HMM_Island")]))
included_annotations <- included_annotations[match(rownames(betas), rownames(included_annotations)),]
res_annotated <- cbind(res, included_annotations)


#save the results to a csv as as
write.csv(res_annotated, "EWAS_GSE66351/Annotated_Results_dataset_test.csv")

#save the results in the standard format defined for a BED file - needed for the DMR calling function
# standard format that requires the first three columns to be: chrom, start and end. 
# This information was downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/
# for the method the 4th column is assumed to contain the p-value for a region
standard_bed_columns <- read.table("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\HAIB.A549.EtOH.Rep.3.bed")
standard_bed_columns <- standard_bed_columns[ ,c(1,2,3,4)]

res_bed <- merge(standard_bed_columns, res, by.x = 4, by.y = "row.names", all.y = TRUE)

res_bed <- data.frame(chrom = res_bed[2], start = res_bed[3], stop = res_bed[4], p_value = res_bed$Diagnosis_P, coeffi = res_bed$Diagnosis_Beta, stan_er = res_bed$Diagnosis_SE, Illumina_ID = res_bed[1])
write.table(res_bed, "EWAS_GSE66351/results_test.bed")
