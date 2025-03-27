import pandas as pd
import numpy as np
import subprocess
import re
import os
import sys

import methylcheck
try:
    from Federated_Differential_Methylation_Analysis.Python_translation import dasen_normalisation
except ModuleNotFoundError:
    sys.path.append("/cosybio/project/vanElferen/FedEWAS")
    from Federated_Differential_Methylation_Analysis.Python_translation import dasen_normalisation

from sklearn.decomposition import PCA 
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
try:
    from Federated_Differential_Methylation_Analysis.Python_translation import EWAS_central_Parallel
except ModuleNotFoundError:
    sys.path.append("/cosybio/project/vanElferen/FedEWAS")
    from Federated_Differential_Methylation_Analysis.Python_translation import EWAS_central_Parallel
try:
    from Federated_Differential_Methylation_Analysis.Evaluations.CreateDesignMatrices import createDesignMatrix66351, createDesignMatrix105109, createDesignMatrix134379
except ModuleNotFoundError:
    sys.path.append("/cosybio/project/vanElferen/FedEWAS")
    from Federated_Differential_Methylation_Analysis.Evaluations.CreateDesignMatrices import createDesignMatrix66351, \
        createDesignMatrix105109, createDesignMatrix134379
import argparse

designMatricesFunctions = {'GSE66351': createDesignMatrix66351,'GSE105109': createDesignMatrix105109, 'GSE134379': createDesignMatrix134379,
                           'GSE66351_half': createDesignMatrix66351,'GSE105109_half': createDesignMatrix105109, 'GSE134379_half': createDesignMatrix134379}

parser = argparse.ArgumentParser(description = "In case of raw data run centralised preprocessing, perform dasen normalisation, run linear model EWAS and save results")
parser.add_argument("input_dir", metavar = "input_dir/", type=str, nargs=1,
    help = "Directory where the *pheno.csv, manifest.csv, /idat are located in case of raw data. Directory where the Filtered_methylated_intensities.csv and Filtered_unmethylated_intensities.csv file are located in case of filtered data.")
parser.add_argument("output_dir", metavar="output_dir/", type = str, nargs=1,
    help = "The basename of the directory where the output directories/files should be created")
parser.add_argument("-i", "--identifier", type=str, nargs=1,
    help = "the identifier used in the output file structure create to save the output of the preprocessing, normalisation and EWAS")
parser.add_argument("--Filtered", "-f", action="store_true",
    help = "Indicate that input files contain already preprocessed and filtered intensities. In case of raw input data leave this out")
parser.add_argument("-scr_dir", "--script_directory", metavar= "script_directory/", type=str, nargs=1,
    help = "Pass the filepath to the directory containing the scirpts needed for preprocessing if they are not in the current working directory")

args = parser.parse_args()
print(args)
input_dir = args.input_dir[0]
output_dir = args.output_dir[0]
if type(args.script_directory[0]) == str:
    script_dir = args.script_directory[0]
else:
    script_dir = sys.path[0] #set script_dir to the directory from where this script was started
    print("The output can be found here:", script_dir)
identifier = args.identifier[0]
if args.Filtered:
    # read in the centrally r-processed filtered methylated and unmethylated intensities
    if '_half' in identifier:
        pheno = pd.read_csv(os.path.join(preprocessing_results_dir, "Reduced_Pheno_Info.csv"), index_col=0)
    else:
        pheno = pd.read_csv(os.path.join(input_dir, (identifier + '_pheno.txt')), index_col=0, sep='\t').T
    unmeth = pd.read_csv(os.path.join(input_dir, "Filtered_Unmethylated.csv"), index_col=0)
    unmeth.astype(np.float64)
    meth = pd.read_csv(os.path.join(input_dir, "Filtered_Methylated.csv"), index_col=0)
    meth.astype(np.float64)
    pre_norm_betas = pd.read_csv(input_dir, "Filtered_Betas.csv", index_col=0)
    pre_norm_betas.astype(np.float64)
    # create the output file structure
    output_dir_QC = os.path.join(output_dir, "QC_Python")
    if not os.path.isdir(output_dir_QC):
        os.makedirs(output_dir_QC)
else:
    try:
        rscript_path = subprocess.check_output("which Rscript", shell=True).decode("utf-8").strip()
        if '_half' in identifier:
            identifier_temp = identifier.split('_')[0]
            preprocessing = subprocess.run([rscript_path, '--vanilla', os.path.join(script_dir, "centralised_preprocessing_half.r"),
            os.path.join(input_dir, "idat"), os.path.join(input_dir, (identifier_temp + "_pheno.txt")), output_dir,
            os.path.join(input_dir, "GPL13534_HumanMethylation450_15017482_v.1.1.csv"), identifier], capture_output=True)
        else:
            preprocessing = subprocess.run(
                [rscript_path, '--vanilla', os.path.join(script_dir, "centralised_preprocessing.r"),
                 os.path.join(input_dir, "idat"), os.path.join(input_dir, (identifier + "_pheno.txt")), output_dir,
                 os.path.join(input_dir, "GPL13534_HumanMethylation450_15017482_v.1.1.csv"), identifier],
                capture_output=True)
        print(preprocessing.stderr)
        print('preprocessing complete')
    except:
        print('There was an error in preprocessing')
        print(preprocessing.stderr)
        sys.exit(1)
    preprocessing_result_dir = os.path.join(output_dir, ("QC_" + identifier))
    if '_half' in identifier:
        pheno = pd.read_csv(os.path.join(preprocessing_result_dir, "Reduced_Pheno_Info.csv"), index_col=0)
    else:
        pheno = pd.read_csv(os.path.join(preprocessing_result_dir, "Pheno_Info.csv"), index_col=0)
    unmeth = pd.read_csv(os.path.join(preprocessing_result_dir, "Filtered_Unmethylated.csv"), index_col=0)
    unmeth.astype(np.float64)
    meth = pd.read_csv(os.path.join(preprocessing_result_dir, "Filtered_Methylated.csv"), index_col=0)
    meth.astype(np.float64)
    pre_norm_betas = pd.read_csv(os.path.join(preprocessing_result_dir, "Filtered_Betas.csv"), index_col=0)
    pre_norm_betas.astype(np.float64)
    # create the output file structure
    output_dir_QC = os.path.join(output_dir, "QC_Python")
    if not os.path.isdir(output_dir_QC):
        os.makedirs(output_dir_QC)


annotation = pd.read_csv(os.path.join(input_dir, "GPL13534_HumanMethylation450_15017482_v.1.1.csv"), skiprows=7, index_col=0, low_memory=False)
annotation_data = annotation.loc[list(set(meth.index.values).intersection(set(annotation.index.values))), :]

# pre-norm betas distribution density plot
pre_norm_BD = methylcheck.beta_density_plot(pre_norm_betas, return_fig=True)
pre_norm_BD.savefig(os.path.join(output_dir, identifier + "_pre-normalisation beta density distribution.jpeg"))

print("Normalising data")
# dasen normalisation
normalised_betas = dasen_normalisation.dasen_normalisation(unmeth, meth, annotation_data.loc[:,"Infinium_Design_Type"])
# save the normalised betas for method comparison
normalised_betas.to_csv(os.path.join(output_dir_QC, (identifier + "_normalised_betas_python.csv")))

# post-norm betas distribution density plot
post_norm_BD = methylcheck.beta_density_plot(normalised_betas, return_fig=True)
post_norm_BD.savefig(os.path.join(output_dir_QC, (identifier + "_post-normalisation beta density distribution.jpeg")))
# PCA to check for batch effects
print("PCA")
n_comp = 10
PCA = PCA(n_components=n_comp)
PCA_out = PCA.fit(normalised_betas)
Principle_components = pd.DataFrame(PCA_out.components_.T, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"], index=pheno.index)
if "Sentrix_ID" in pheno.columns:
    sentrixCol = "Sentrix_ID"
elif "sentrix_id" in pheno.columns:
    sentrixCol = "sentrix_id"
else:
    print('there is no sentrix_id information in the pheno information, generating it from sample barcodes')
    pheno['Sentrix_ID'] = [i.split('_')[1] for i in pheno.Sample_ID]
    sentrixCol = 'Sentrix_ID'
    print('saving pheno information with Sentrix_ID column')
    if '_half' in identifier:
        pheno.to_csv(os.path.join(preprocessing_result_dir,"Reduced_Pheno_Info.csv"))
    else:
        pheno.to_csv(os.path.join(preprocessing_result_dir, "Pheno_Info.csv"))
data_to_plot = pd.concat([Principle_components, pheno.loc[:, sentrixCol]], axis=1)
# exploratory grid plot to visualise the correlation between the PCs and Sentrix_ID
g = sns.PairGrid(data_to_plot, hue= sentrixCol)
g.map(sns.scatterplot)
try:
    g.savefig(os.path.join(output_dir_QC, (identifier + "_Exploratory PCA grid coloured by Sentrix_ID.jpeg")))
except:
    print("Issue with saving figure")
# calculate the correlation between the PC clusters and the Sentrix_ID variable
pca_iterator = np.arange(0, n_comp)
pca_cor = []
codes_sentrix, options_sentrix = pd.factorize(pheno.loc[:,sentrixCol])
for id, x in np.ndenumerate(pca_iterator):
    cor, p_value = pearsonr(PCA_out.components_[x,], codes_sentrix)
    pca_cor.append(float(cor)) # turn into list append
pc_output = pd.DataFrame({"PC": pca_iterator, "Correlation": pca_cor})
pc_output.sort_values(by = "Correlation", key = pd.Series.abs, ascending = False, inplace=True) #sort the correlation values in descending order based on their absolute value
first = abs(pc_output.iloc[0]) #get the two highest correlated PCs
second = abs(pc_output.iloc[1])
print("The top correlated principle components with Sentrix_ID:", first.name + 1, second.name + 1)

# visualise the two PCs with the strongest correlation with Sentrix_ID
data_to_plot = pd.DataFrame({"PC4":PCA_out.components_[int(first.name)], "PC5":PCA_out.components_[int(second.name)], "Sentrix_ID":pheno.loc[:, sentrixCol]},
     index=pheno.index.values)
detail_PCA = sns.lmplot(x="PC4", y = "PC5", hue="Sentrix_ID", data = data_to_plot, fit_reg=False)
try:
    detail_PCA.savefig(os.path.join(output_dir_QC, (identifier + "_Highest correlating PCs with Sentrix_ID coloured by Sentrix_ID.jpeg")))
except:
    print("Issue with saving figure")

# EWAS
# check if the design matrix exists
if not os.path.exists(os.path.join(input_dir, "Small_EWAS_design.csv")):

    if '_half' in identifier:
        designMatricesFunctions[identifier](pheno_df_path=os.path.join(preprocessing_result_dir,"Reduced_Pheno_Info.csv"),
                            small=True, federated=False, per_region=False,
                            output_path=input_dir)
    else:
        designMatricesFunctions[identifier](pheno_df_path=os.path.join(preprocessing_result_dir,"Pheno_Info.csv"),
                                small=True, federated=False, per_region=False,
                                output_path=input_dir)

design_matrix = pd.read_csv(os.path.join(input_dir, "Small_EWAS_design.csv"), index_col=0)

print("EWAS")
results_diagnosis, results_ewas = EWAS_central_Parallel.EWAS_central(design_matrix, normalised_betas, identifier)

# create an output table with the top (genomewide) significant probes and their associated gene with the metrics
gene_annotations = annotation.loc[:, "UCSC_RefGene_Name"]
gene_annotations = gene_annotations.loc[set(results_ewas.index.values).intersection(set(annotation_data.index.values))]
final_results_table = pd.DataFrame({"Associated Gene":gene_annotations, "Methylation Change":results_ewas.loc[:,("Coefficient", "AD")], "Corrected P-value":results_ewas.loc[:,("Corrected P-value", "AD")]})
final_results_table.sort_values(by = ["Corrected P-value"], inplace = True)
# save the EWAS results to the output directory
final_results_table.to_csv(os.path.join(output_dir_QC, (identifier + "_results_diagnosis_regression_python.csv")))
results_ewas.to_csv(os.path.join(output_dir_QC, (identifier + "_full_results_regression_python.csv")))
# save the SSE to a file for the calculation of the eBayes results
with open(os.path.join(output_dir_QC, (identifier + "_SSE_list.csv")), "w") as outfile:
    outfile.writelines('\n'.join([str(i) for i in SSE]))
# try the ebayes calculation with the regression output
regressioResultsCalculator = eBayesLocal.eBayesLocal(results_ewas,
                                         SSE,
                                         design_matrix.shape[0])
regressioResultsCalculator.eBayes()
regressioResultsCalculator.table.to_csv(os.path.join(output_dir_QC, (identifier + "_eBayesTopTableResult.csv")))

sys.exit(0)


    
