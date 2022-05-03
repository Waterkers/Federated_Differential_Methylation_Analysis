input <- commandArgs(trailingOnly = TRUE)

msetEPIC.pf <- input[1]
QCmetrics <- input[2]
manifest_path <- input[3]
output_dir <- input[4]


install.packages("https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/RefFreeEWAS_2.2.tar.gz", repos = NULL, type = "source") # not available on CRAN anymore so needs to be downloaded manually from the archive
# add the filepath to where the archived version of RefFreeEWAS was downloaded above
install.packages("quadprog", repos = "https://mirror.lyrahosting.com/CRAN/") # one of the RefFreeEWAS dependencies that I hadn't installed yet
library(RefFreeEWAS)
need <- c("wateRmelon", "methylumi")
if (!require(need, quietly = TRUE))
  BiocManager::install(need)
library(wateRmelon, methylumi)
install.packages("tidyverse", repos = "https://mirror.lyrahosting.com/CRAN/")
library(tidyverse)

cell_decomp_RefFreeEWAS <- function(manifest_path, msetEPIC.pf, QCmetrics, output_dir){

# load in the methylSet object
	load(msetEPIC.pf)
	pheno <- QCmetrics
	annotation_data <- read.csv(manifest_path , skip = 7, header = TRUE)
	retained_annotation <- annotation_data[annotation_data$IlmnID %in% rownames(methylumi::betas(msetEPIC.pf)), ]
	fData(msetEPIC.pf)["CHR"] <- retained_annotation["CHR"]

	

##### Cell type estimation ##############
# start with removing the X-chromosome probes from the dataset if they are there
# because they can interfere with the cell-type decomposition.
	temp_data = betas(msetEPIC.pf[fData(msetEPIC.pf)$CHR!="X", ]) #something funky going on here


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
# plot(elbowdif)
# abline(h=quantile(elbowdif,0.98),col="red")
# points(elbowdif[elbowdif>quantile(elbowdif,0.98)],pch=19,col="blue")

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

	pdf(file.path(output_dir, "Celltype_estimates_barplot.pdf"))
	ggplot(DF, aes(x = SampleID, y = Proportion, fill = Celltypes))+
		ggtitle("Celltype composition estimate per sample")  + 
		coord_flip() +
		geom_bar(stat = "identity")
	dev.off()

#### save the final pre-processed betas with minimal and full phenotype information
	Betas <- betas(msetEPIC.pf)

	temp_Pheno <- QCmetrics[match(colnames(Betas), QCmetrics$Sample_ID), ]
# remove unnecessary text from the phenotype dataframe cells
	temp_Pheno <- lapply(temp_Pheno, sub, pattern = "^[^:]*:", replacement = "")

	Cell_Types <- CT[match(colnames(Betas), rownames(CT)),]

	Small_Pheno <- data.frame(Sample_ID = temp_Pheno$Sample_ID, Diagnosis = temp_Pheno$Diagnosis, Sex = temp_Pheno$Sex,
                          Age = temp_Pheno$Age, Cell_Type = Cell_Types)

# Create the full phenotype file
	Full_Pheno <- data.frame(temp_Pheno, Sample_ID = QCmetrics$Sample_ID, Cell_Type = Cell_Types)

# save everything
	save(Betas, Small_Pheno, file = file.path(output_dir, "Post_RefFreeEWAS.RData"))

	save(Full_Pheno, file = file.path(output_dir, "Full_Phenotype_Information.RData"))

#save the information stored in the data2 object (the MethylumiSet object thing) into seperate dataframes
	
# betas
	write.csv(Betas, file.path(output_dir, "Final_preprocessed_betas.csv"))

# full phenotype information
	write.csv(Full_Pheno, file.path(output_dir, "Full_pheno_information.csv"))

	print("Preprocessed data, short and full phenotype information saved")
}

cell_decomp_RefFreeEWAS(manifest_path, msetEPIC.pf, QCmetrics, output_dir)