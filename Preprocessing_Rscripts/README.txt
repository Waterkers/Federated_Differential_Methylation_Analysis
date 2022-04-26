This folder contains the r-scripts detailing the data preprocessing pipeline supplied by <group> at Exeter University.
The order of preprocessing steps has been modified based on own insight:
	- The removal of SNP and cross-hybridising probes has been moved up before normalisation of the data. 
	- Cell type decomposition code was added to the complete preprocessing script which was not supplied originally.
		This is based on an archieved version of RefFreeEWAS downloaded from:
		
The scripts in this folder can be used for:
	centralised_preprocessing.r contains the steps of the preprocessing pipeline that are to be preformed locally before
	loading the data into the federated implementation of the dasen normalisation function.
		Loading idat files
		Checking for absolute outliers in terms of (un)methylated signal intensities
		Checking bisulfite conversion efficiency for the samples
		Checking for correct labeling of biological sex of the samples
		Checking for correlation between SNPs - genetically identical samples
		P-filter (wateRmelon) remove probes with low signal intensity and samples with low quality data
		Remove cross-hybridising and known SNP associated probes - based on the McCartney et al., 2016 paper
	The files saved by this script are:
		Raw_betas.csv
		Raw_methylated_intensities.csv
		Raw_unmethylated_intensities.csv
		Preprocessed_betas.csv
		Preprocessed_methylated_intensities.csv
		Preprocessed_unmethylated_intensities.csv
		pre_norm_pheno_information.csv
	The .RData files saved by this script:
		methylumiSet_PhenoDatafrma.RData - MethylumiSet object with raw data and the accompanying phenotype file
		preprocessed_MethyLumiSet.RData - MethylumiSet object with preprocessed data and the accompanying phenotype file
	
	Loading_idats_code_saveOutput_python_shell.R only reads in idat files and saves the output in .csv files as well as an .RData file
	The output .RData file also contains the phenotype information with the Sample_ID column created based on the methylumi GetBarcodes function.
	The .csv files saved by this script are:
		methylated_intensities.csv
		unmethylated_intensities.csv
		betas.csv
		phenotype.csv
	The .RData files saved by this script
		methylumiSet_PhenoDatafrma.RData - MethylumiSet object with raw data and the accompanying phenotype file
	
	