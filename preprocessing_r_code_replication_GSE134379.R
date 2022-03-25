##### Start with installing the required packages ########
need <- c("wateRmelon", "methylumi", "ChAMP")
if (!require(need, quietly = TRUE))
  BiocManager::install(need)
library(wateRmelon, methylumi)
library(ChAMP) # for some reason library only loads this package when it is called on its own

#if (!require("RefFreeEWAS", quietly = TRUE))
install.packages("C:\\Users\\Silke\\OneDrive\\Documenten\\R\\win-library\\4.1\\RefFreeEWAS_2.2.tar.gz", repos = NULL) # not available on CRAN anymore so needs to be downloaded manually from the archive
# for the workflow probably needs full translation
install.packages("quadprog") # one of the RefFreeEWAS dependencies that I hadn't installed yet
library(RefFreeEWAS)
##### source the Exeter functions needed for the pipeline
lapply(list.files("F:\\Msc Systems Biology\\MSB5000_Master_Thesis\\ExeterEWASPipeline-master\\R",pattern = "\\.r$",full.names = T),function(x){source(x)})
## Set the working directory
setwd("F:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\R_Pipeline") # set working directory
## create a folder for the QC output
if(!dir.exists("QC_GSE134379")){
  dir.create("QC_GSE134379")
}
if(!dir.exists("QC_GSE134379/Plots")){
  dir.create("QC_GSE134379/Plots")
}
#### import the data using methylumi -> create a MethyLumiSet object

# loading in actual data - GSE134379
pheno1 <- read.table("F:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE134379_RAW\\GSE134379_pheno.txt")
pheno1 <- t(pheno1)#transpose the imported tabel to the sample characteristics/ids etc are columns and the samples are rows
pheno1 <- as.data.frame(pheno1)
colnames(pheno1)<- pheno1[1,]
pheno1 <- pheno1[2:191,]
Sample_ID <- getBarcodes("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE134379_RAW\\idat")
pheno1 <- cbind(pheno1, Sample_ID)

# set up the system for parallel processing to make it possible to deal with the dataset
install.packages("parallel")
library(parallel)
num_cores <- detectCores()

#data1 <- methylumi::methylumIDAT(barcodes = Sample_ID, pdat = pheno1, idatPath = "E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE134379_RAW\\idat", parallel = TRUE, mc.cores = num_cores)
data2 <- wateRmelon::readEPIC(barcodes = Sample_ID, pdat = pheno1, idatPath = "E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE134379_RAW\\idat", parallel = TRUE, mc.cores = num_cores)
#fData(data1) <- read.csv("F:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE66351_RAW\\GPL13534_HumanMethylation450_15017482_v.1.1_edit.csv", header = TRUE)
## next it is necessary to rename the phenotype and data object to the names that are used in the pipeline
pheno <- pheno1
msetEPIC <- data2 #or data1 depending on the read in method used

########## Start with the QC pipeline from Exeter ##################
### checking methylated and unmethylated intensities #############
### extract sample intensities 
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
intens.Thres<-2000 

# make PDF histogram of sample intensities
pdf("QC/Plots/Sample_Intensity_histogram.pdf")
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


## plot each plate as a boxplot - this does not work for the sample data since there is only a small set
# of samples provided
pdf("QC/Plots/Sample_Intensity_ByPlate_boxplot.pdf")
par(mfrow = c(1,2))
par(mar = c(8, 4, 1, 1))
nCol<-length(unique(pheno$plate))## assumes there is a column called Plate in your phenotype file
boxplot(M.median ~ pheno$plate, ylab = "Median M intensity", xlab = "Plate", las = 2, col = rainbow(nCol)) 
boxplot(U.median ~ pheno$plate, ylab = "Median U intensity", xlab = "Plate", las = 2, col = rainbow(nCol))
dev.off()
# 
# ## alternatively colour points in original scatterplot by Plate
# pdf("QC/Plots/Sample_Intensity_ByPlate_histogram.pdf")
# nCol<-length(unique(pheno$plate))## assumes there is a column called Plate in your phenotype file
# plot(M.median, U.median, pch = 16, xlab = "Median M intensity", ylab = "Median U intensity", col = rainbow(nCol)[factor(pheno$plate)])
# abline(v = intens.Thres, col = "red")
# abline(h = intens.Thres, col = "red") 
# legend("topright", levels(factor(pheno$plate)), col = rainbow(nCol), pch = 16)
# dev.off()

##### Checking the Bisulfite conversion ################
bs<-wateRmelon::bscon(msetEPIC) # - doesn't work because of this error "Error in bsI.green[1:2, ] : subscript out of bounds"

pdf("QC/Plots/Bisulphite_Conversion.pdf")
hist(bs, xlab = "Median % BS conversion", main = "")
abline(v = 80, col = "red")
dev.off()
QCmetrics<-cbind(QCmetrics, bs)

##### Checking for biological sex vs. measured sex ############
betas <- methylumi::betas(msetEPIC)
pheno<-pheno[match(colnames(betas), pheno$Sample_ID),]

pdf("QC/Plots/Gender_Cluster.pdf")
predSex1<-findGenderPC(betas, pheno$Sample_sex, npcs = 20) # the default setting of npcs = 20 gave an error so
# I reduced the number of princicple components for the example data set -> remember to put it back at 20
# when using real data
predSex2<-clusterGender(betas, pheno$Sample_sex)
dev.off()

# Confirm findings
PCA = prcomp(betas[complete.cases(betas),])
plot(PCA$rotation[,1],PCA$rotation[,3],col=ifelse(pheno$Sample_sex=="Sex: M","blue","magenta"), main= "PCA of betas, colored by phenotype trait sex",pch=19)
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

pdf("QC/Plots/SNP_Correlations.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()

# Check which samples match to their best match 
pdf("QC/Plots/SNP_Samples.pdf", width = 15, height = 8)
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
  if(!id_match){warning(paste0("mismatch with sample ",xname, " Index: ",i))}
  
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
pdf("QC/Plots/Betas_raw_boxplot.pdf")
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

###### Normalisation ################
msetEPIC.pf <- dasen(msetEPIC.pf)
# adding feature/probe annotation information after normalisation because otherwise
# the dasen function becomes fussy and won't work
annotation_data <- read.csv("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\GSE66351_RAW\\GPL13534_HumanMethylation450_15017482_v.1.1_edit.csv", header = TRUE)
retained_annotation <- annotation_data[annotation_data$Probe_ID %in% rownames(betas(msetEPIC.pf)), ]
fData(msetEPIC.pf) <- retained_annotation

pdf("QC/Plots/Betas_normalized_boxplot.pdf")
boxplot(betas(msetEPIC.pf),main="Betas normalized")
dev.off()

##### Cell type estimation ##############
# start with removing the X-chromosome probes from the dataset if they are there
# because they can interfere with the cell-type decomposition.
temp_data = betas(msetEPIC.pf[fData(msetEPIC.pf)$Chromosome_36!="X", ]) #something funky going on here


# next transpose the data so the samples become the rows and run signular value decomposition
# the top x most informative singular value decomposition 'components' are used to select the most
# informative probes to be used for the cell-type decomposition. 
sv = svd(scale(t(temp_data)))    # PCA ; samples should become row
plot(sv$d^2,type="b")        # Scree plot - shows the most informative data decompositions based on the input data

# from here the most informative decompositions/values/whatever they are called need to be selected. The top 98th
# percentile of values is used here
# values for elbow to detect
elbowdif = sv$d^2

# estimate better amount of PCs (based on above 98th percentile)
s_maxPC=10
s_minPC=1
s_maxPC = min(s_maxPC,max(s_minPC,sum(elbowdif>quantile(elbowdif,0.98)))) # calculate the 98th percentile based on the
# svd output calculated based on the input data and store the number of values that are above this cut-off in the term
# s_maxPC

# Automated elbow cutoff - visualisation of the values selected (above the cut-off)
plot(elbowdif)
abline(h=quantile(elbowdif,0.98),col="red")
points(elbowdif[elbowdif>quantile(elbowdif,0.98)],pch=19,col="blue")

# select the right singular vectors corresponding to the selected top values
pcSelect =  sv$v[,1:s_maxPC]  
pcSelect = as.data.frame(pcSelect)
# from here select the probes with the most informative factor-loadings - using mahalanobis distance
mh = mahalanobis(pcSelect, apply(pcSelect,2,mean), cov(pcSelect)) # calculate mahalanobis distance for all the 
# factor-loadings
hist(mh) #plot a histogram of these

#select the top 5000 (as initial estimate) with the lowest mh distance -> these are most informative
s_select_CpGs = 10000 # first 5000 OK
plot(mh[order(mh, decreasing=TRUE)])
abline(v=s_select_CpGs,col="red",lty=2) # plot the mh distances in decreasing order together with the suggested
# cut off for informative CpGs to see if all of the factor loadings with a mh distance < 1 are included
# now select the most important CpGs based on the cut-off defined above
cpgSelect = order(mh, decreasing=TRUE)[1:s_select_CpGs]
Ysel = temp_data[cpgSelect,]
rm(temp_data) # clear space in memory
# as a final step before moving on the cell type decomposition is to double check that there are no sex effects
# present in the data anymore
plot((prcomp(t(Ysel))$x),col=ifelse(pheno$Sample_Sex=="Sex: M",1,2))

# Maximum number of estimatable celltypes set to 5 (this is default and supported by evidence see https://www.nature.com/articles/s41598-018-25311-0 )
s_maxCelltypes  = 5

# Get PCs (without standardization)
svSel2 = svd(Ysel)

# Initial deconvolution (note, this could take awhile)
cellmixArray  <- RefFreeCellMixArrayWithCustomStart(Ysel, 
                                                    mu.start = svSel2$u,  # Initial methylome matrix 
                                                    Klist=3:s_maxCelltypes          # List of K values to try (# constituent cell types)
)

rm(svSel2) 
gc()

#-----------------------------------------------------------------------------------------------------#
#                           Bootstrap
#-----------------------------------------------------------------------------------------------------#
# Do the bootstrap for selecting the K parameter (# assumed cell types)
cellmixArrayBoot <- RefFreeCellMixArrayDevianceBoots(
  cellmixArray,            # Array object
  Y=Ysel,                  # Data on which array was based
  R=100,                    # defaults 5; 
  bootstrapIterations=5)    # defaults 5; 

#-----------------------------------------------------------------------------------------------------#
#                           Select K
#-----------------------------------------------------------------------------------------------------#
# Show mean deviance per K, per window
wnsrMeanDev <-apply(cellmixArrayBoot[-1,], 2, mean, trim=0.25)
#wnsrMeanDev

# Choose K based on minimum deviance
Kchoose <- as.numeric(which.min(wnsrMeanDev))
#Kchoose # 1

#-----------------------------------------------------------------------------------------------------#
#                           Omega
#-----------------------------------------------------------------------------------------------------#
# Chosen Omega
Omega <- cellmixArray[[ Kchoose ]]$Omega # The cell mixture matrix
MuSmall <- cellmixArray[[ Kchoose ]]$Mu

# Celltype estimations (Omega matrix)
Omega = data.frame(Omega)
colnames(Omega)=paste0("CT",1:(dim(Omega)[2]))

# rename to CT
CT= Omega
rm(Omega)

# make image
temp_CT = CT
temp_CT$SampleID=rownames(CT)

# format data
DF <- pivot_longer(temp_CT,cols = 1:dim(CT)[2],names_to = "Celltypes",values_to = "Proportion") %>% group_by(SampleID,Celltypes) #%>%  mutate(name = fct_reorder(name, value)) 

pdf("QC/Plots/Celltype_estimates_barplot.pdf")
ggplot(DF, aes(x = SampleID, y = Proportion, fill = Celltypes))+
  ggtitle("Celltype composition estimate per sample")  + 
  coord_flip() +
  geom_bar(stat = "identity")
dev.off()

#### save the final pre-processed betas with minimal and full phenotype information
Betas <- betas(msetEPIC.pf)

temp_Pheno <- QCmetrics[match(colnames(Betas), QCmetrics$Sample_ID)]
# remove unnecessary text from the phenotype dataframe cells
temp_Pheno <- lapply(temp_Pheno, sub, pattern = "^[^:]*:", replacement = "")

Cell_Types <- CT[match(colnames(Betas), rownames(CT)),]

Small_Pheno <- data.frame(Sample_ID = QCmetrics$Sample_ID, Diagnosis = temp_Pheno$Sample_diagnosis, Sex = temp_Pheno$Sample_sex,
                          Age = temp_Pheno$Sample_age, Cell_Type = Cell_Types)
# create a column with the full sentrix ID because it seems handy
Small_Pheno$Sentrix_ID <- str_c(temp_Pheno$Sample_sentrix_id, "_", temp_Pheno$Sample_sentrix_position)

# Create the full phenotype file
Full_Pheno <- data.frame(temp_Pheno, Sample_ID = QCmetrics$Sample_ID, Cell_Type = Cell_Types)

# save everything
save(Betas, Small_Pheno, file = "QC\\GSE66351_first20_QCandDasen.RData")
save(Full_Pheno, file = "QC\\GSE66351_Full_Phenotype_Information.RData")

