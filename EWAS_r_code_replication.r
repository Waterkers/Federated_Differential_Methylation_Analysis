############################################################################
#                         EWAS                                             #
###########################################################################
# load the exeter functions
lapply(list.files("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\ExeterEWASPipeline-master\\R",pattern = "\\.r$",full.names = T),function(x){source(x)})
# setwd("") # set the working directory
load("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\R_Pipeline\\QC\\GSE66351_first20_QCandDasen.RData")
## create a folder to save EWAS output
if(!dir.exists("EWAS")){
  dir.create("EWAS")
}
if(!dir.exists("EWAS/Plots")){
  dir.create("EWAS/Plots")
}

#### Get rid of cross-hybridising probes, based on McCartney et al., 2016 ###
crosshyb <- read.table("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Federated_Differential_Methylation_Analysis\\Cross_hybridising_CpGTargetting_Probes_McCartneyetal2016.txt")

betas<-Betas[!(rownames(Betas) %in% crosshyb[,1]), ]
betas<-filterSNPprobes(betas, population = "EUR", maf = 0.05) ## filters common probes based on allele frequency in european populations.
betas<-betas[-grep("rs", rownames(betas)),] ## remove SNP probes