Complete r analysis workflow

The complete 450k methylation data analysis workflow used in this project is included in this folder.
It consists of two scripts:
	complete_preprocessing.r
	EWAS.r

Each script starts with a section that indicates which variables, filepaths etc. need adapting by the user 
before use.

complete_preprocessing.r contains the following preprocessing steps:
		Reading in the data using the readEPIC function from the wateRmelon package
		Investigating the measured methylated and unmethylated intensities for each probe grouped by sample
			and grouped by plate if this information is present in the phenotype information
		Checking the bisulfite conversion efficiency for the samples
		Checking if the recorded sex in the phenotype information matches the measured sex based on the data
		Checking for correlation between known SNP related probes to detect any genetically identical samples
		Filtering probes and samples based on signal quality
		Removing cross-hybridising probes and probes that are known to be part of genes containing common SNPs 
		in the population under investigation (the probe lists can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830)
		Normalising the data kept after the two filtering steps using dasen normalisation (wateRmelon package)
		Performing cell type decomposition using RefFreeEWAS version 2.2 (archived version available here: https://cran.r-project.org/src/contrib/Archive/RefFreeEWAS/)
	This script takes the following as input:
		a phenotype information file 
		.idat files (number of red/green file pairs == number of samples in the phenotype information file)
		probe annotation file
	This script produces the following output when run:
		.csv files:
			Raw_Betas.csv
			Raw_Methylated.csv
			Raw_Unmethylated.csv
			Filtered_Betas.csv
			Filtered_Methylated.csv
			Filtered_Unmethylated.csv
			Normalised_Betas.csv
			Preprocessed_CTD_Betas.csv
			Full_Pheno_Info.csv
			Reduced_Pheno_Info.csv
		.RData files:
			FilteredMethylumisetQCmetrics.RData
			Normalised.RData
			Preprocessed_CTD.RData - includes the information saved in Reduced_Pheno_Info.csv
			Full_Phenotype_Information.RData

EWAS.r contains a linear model that calculates the effect the regression coefficients for each probe and independent variable combination to provide
an indication of the extend to which the probes can account for changes in the independent variables. 
	The script takes the following input:
		Matrix or dataframe of beta values
		Dataframe of phenotype information to be used in the linear model. Columns used in the script:
			Diagnosis
			Sex
			Age
			Sentrix_ID
			CellType.CTx - output of the RefFreeEWAS cell type decomposition
	
	The script produces the following output:
		Results_dataset.csv - contains the coefficient, standard error and p-values for each probe for the diagnosis (first independent variable in the model)
		Annotated_Results_dataset.csv - annotated version of the results dataset
		results.bed - the results from the linear model presented in a .bed format
	The independent variables for which the linear model results are saved should be modified based on the research question
		