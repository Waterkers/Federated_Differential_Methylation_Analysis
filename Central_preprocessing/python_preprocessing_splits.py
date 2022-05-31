from fnmatch import fnmatch
import pandas as pd
import numpy as np
import subprocess
import re
import os
import sys
import fnmatch

import methylcheck
import dasen_normalisation

from sklearn.decomposition import PCA 
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import EWAS_central

import argparse

parser = argparse.ArgumentParser(description = "In case of raw data run centralised preprocessing, perform dasen normalisation, run linear model EWAS and save results")
parser.add_argument("input_dir", metavar = "input_dir/", type=str, nargs=1,
    help = "Directory where the *pheno.csv, manifest.csv, /idat are located in case of raw data. Directory where the Filtered_methylated_intensities.csv and Filtered_unmethylated_intensities.csv file are located in case of filtered data.")
parser.add_argument("output_dir", metavar="output_dir/", type = str, nargs=1,
    help = "The basename of the directory where the output directories/files should be created")
parser.add_argument("-i", "--identifier", type=str, nargs=1,
    help = "the identifier used in the output file structure create to save the output of the preprocessing, normalisation and EWAS")



args = parser.parse_args()
print(args)
input_dir = args.input_dir[0]
output_dir = args.output_dir[0]

identifier = args.identifier[0]

output_dir_QC = os.path.join(output_dir, "QC_Python_splits")
if not os.path.isdir(output_dir_QC):
    os.makedirs(output_dir_QC)

files = [file for file in os.listdir(input_dir) if fnmatch.fnmatch(file, "Split*")]
file_split = []
for file in files :
    file_split.append(re.search("(Split_[1-9])", file).group(1))
unique_splits = set(file_split)
for split in unique_splits:
    meth = pd.read_csv(os.path.join(input_dir, split + "_methylated.csv"), index_col=0)
    unmeth = pd.read_csv(os.path.join(input_dir, split + "_unmethylated.csv"), index_col=0)
    pre_norm_betas = pd.read_csv(os.path.join(input_dir, split + "_betas.csv"), index_col=0)
    pheno = pd.read_csv(os.path.join(input_dir, split + "_pheno.csv"), index_col=0)
    design_matrix = pd.read_csv(os.path.join(input_dir, split + "_design.csv"), index_col=0)
    

    annotation = pd.read_csv(os.path.join(input_dir, "GPL13534_HumanMethylation450_15017482_v.1.1.csv"), skiprows=7, index_col=0, low_memory=False)
    annotation_data = annotation.loc[set(meth.index.values).intersection(set(annotation.index.values)), :]

    # pre-norm betas distribution density plot
    pre_norm_BD = methylcheck.beta_density_plot(pre_norm_betas, return_fig=True)
    pre_norm_BD.savefig(os.path.join(output_dir_QC, (identifier + split + "_pre-normalisation beta density distribution.jpeg")))

    print("Normalising data")
    # dasen normalisation
    normalised_betas = dasen_normalisation.dasen_normalisation(unmeth, meth, annotation_data.loc[:,"Infinium_Design_Type"])
    # save the normalised betas for method comparison
    normalised_betas.to_csv(os.path.join(output_dir_QC, (identifier + split + "_normalised_betas_python.csv")))

    # post-norm betas distribution density plot
    post_norm_BD = methylcheck.beta_density_plot(normalised_betas, return_fig=True)
    post_norm_BD.savefig(os.path.join(output_dir_QC, (identifier + split + "_post-normalisation beta density distribution.jpeg")))
    # PCA to check for batch effects
    print("PCA")
    n_comp = 10
    PCA = PCA(n_components=n_comp)
    PCA_out = PCA.fit(normalised_betas)
    Principle_components = pd.DataFrame(PCA_out.components_.T, columns=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"], index=pheno.index)
    data_to_plot = pd.concat([Principle_components, design_matrix.loc[:, "Sentrix_ID"]], axis = 1)
    # exploratory grid plot to visualise the correlation between the PCs and Sentrix_ID
    g = sns.PairGrid(data_to_plot, hue= "Sentrix_ID")
    g.map(sns.scatterplot)
    try:
        g.savefig(os.path.join(output_dir_QC, (identifier + split + "_Exploratory PCA grid coloured by Sentrix_ID.jpeg")))
    except:
        print("Issue with saving figure")
    # calculate the correlation between the PC clusters and the Sentrix_ID variable
    pca_iterator = np.arange(0, n_comp)
    pca_cor = []
    codes_sentrix, options_sentrix = pd.factorize(pheno.loc[:,"Sentrix_ID"])
    for id, x in np.ndenumerate(pca_iterator):
        cor, p_value = pearsonr(PCA_out.components_[x,], codes_sentrix)
        pca_cor.append(float(cor)) # turn into list append
    pc_output = pd.DataFrame({"PC": pca_iterator, "Correlation": pca_cor})
    pc_output.sort_values(by = "Correlation", key = pd.Series.abs, ascending = False, inplace=True) #sort the correlation values in descending order based on their absolute value
    first = abs(pc_output.iloc[0]) #get the two highest correlated PCs
    second = abs(pc_output.iloc[1])
    print("The top correlated principle components with Sentrix_ID:", first.name + 1, second.name + 1)

    # visualise the two PCs with the strongest correlation with Sentrix_ID
    data_to_plot = pd.DataFrame({"PC4":PCA_out.components_[int(first.name)], "PC5":PCA_out.components_[int(second.name)], "Sentrix_ID":design_matrix.loc[:, "Sentrix_ID"]},
     index=pheno.index.values)
    detail_PCA = sns.lmplot(x="PC4", y = "PC5", hue="Sentrix_ID", data = data_to_plot, fit_reg=False)
    try:
        detail_PCA.savefig(os.path.join(output_dir_QC, (identifier + split + "_Highest correlating PCs with Sentrix_ID coloured by Sentrix_ID.jpeg")))
    except:
        print("Issue with saving figure")

    # EWAS
    
    print("EWAS")
    results_diagnosis, results_ewas = EWAS_central.EWAS_central(design_matrix, normalised_betas)

    # create an output table with the top (genomewide) significant probes and their associated gene with the metrics
    gene_annotations = annotation.loc[:, "UCSC_RefGene_Name"]
    gene_annotations = gene_annotations.loc[set(results_ewas.index.values).intersection(set(annotation_data.index.values))]
    final_results_table = pd.DataFrame({"Associated Gene":gene_annotations, "Methylation Change":results_ewas.loc[:,("Coefficient", "Diagnosis")], "Corrected P-value":results_ewas.loc[:,("Corrected P-value", "Diagnosis")]})
    final_results_table.sort_values(by = ["Corrected P-value"], inplace = True)
    # save the EWAS results to the output directory
    final_results_table.to_csv(os.path.join(output_dir_QC, (identifier + split + "_results_diagnosis_regression_python.csv")))
    results_ewas.to_csv(os.path.join(output_dir_QC, (identifier + split + "_full_results_regression_python.csv")))

sys.exit(0)


    