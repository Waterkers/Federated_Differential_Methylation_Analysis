#%%
import numpy as np
import os
import pandas as pd

from Federated_Differential_Methylation_Analysis.Central_preprocessing.python_preprocessing_splits import output_dir
from server import Server 
import argparse
from client import Client
#%%
#import pyximport
#pyximport.install(setup_args={"script_args" : ["--verbose"]})
#from linbinR import fast_linbin
#%%
args = None
parser = argparse.ArgumentParser(description='Federated Differential Methylation Analysis')
parser.add_argument('split_directory','-s')
parser.add_argument('output_dir','-o')
parser.add_argument('probe_annotation_path', '-p')
parser.add_argument('split_type', '-t')
args = parser.parse_args()
#%%
if not args:
    split_dir = "/home/silke/Documents/Fed_EWAS/Data/QC_GSE134379_half/GSE134379_splits"
    output = "/home/silke/Documents/Fed_EWAS/Data/QC_GSE134379_half/GSE134379_Fed"
    probeAnnotationPath = "/home/silke/Documents/Fed_EWAS/Data/GSE134379_RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv"
    split_type = 'balanced'
else:
    split_dir = args.split_directory
    output = args.output_dir
    probeAnnotationPath = args.probe_annotation_path
    split_type = args.split_type
#%%
# check if output directory exists, if not make it
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

#%% md
# ## Initialising the clients
#%%
# create client
lab_a = Client("Lab_A", os.path.join(split_dir, "Split_1_design.csv"), os.path.join(split_dir, "Split_1_methylated.csv"), os.path.join(split_dir, "Split_1_unmethylated.csv"), probeAnnotationPath)
lab_b = Client("Lab_B", os.path.join(split_dir, "Split_2_design.csv"), os.path.join(split_dir, "Split_2_methylated.csv"), os.path.join(split_dir, "Split_2_unmethylated.csv"), probeAnnotationPath)
lab_c = Client("Lab_C", os.path.join(split_dir, "Split_3_design.csv"), os.path.join(split_dir, "Split_3_methylated.csv"), os.path.join(split_dir, "Split_3_unmethylated.csv"), probeAnnotationPath)
#%% md
# ## Initialising the server
#%%
serv = Server(["AD", "CTRL"], ["Age", "Sex", "Sentrix_ID"])
global_conditions = serv.return_global_conditions()
#%% md
# ## Joining clients to the server
#%%
# join the clients
serv.get_clients(lab_a.cohort_name, lab_a.probes, lab_a.designmatrix.index)
serv.get_clients(lab_b.cohort_name, lab_b.probes, lab_b.designmatrix.index)
serv.get_clients(lab_c.cohort_name, lab_c.probes, lab_c.designmatrix.index)
#%%
global_probes = serv.find_global_probes()
cohort_effect = serv.find_cohort_effects()
#%%
#check client input
lab_a.input_validation(global_conditions, global_probes)
lab_b.input_validation(global_conditions, global_probes)
lab_c.input_validation(global_conditions, global_probes)
#%%
lab_a.cohort_effects(serv.cohort_effects)
lab_b.cohort_effects(serv.cohort_effects)
lab_c.cohort_effects(serv.cohort_effects)
#%%
if "Sentrix_ID" in global_conditions:
    lab_a.find_unique_SentrixIDS()
    lab_b.find_unique_SentrixIDS()
    lab_c.find_unique_SentrixIDS()
    global_sentrix = serv.return_global_SentrixID(lab_a.unique_SentrixIDS,
                                lab_b.unique_SentrixIDS,
                                lab_c.unique_SentrixIDS)
    lab_a.SentrixID_effects(global_sentrix)
    lab_b.SentrixID_effects(global_sentrix)
    lab_c.SentrixID_effects(global_sentrix)
    


#%% md
# ## Dasen normalisation
#%% md
# Client side
#%%
lab_a.intensity_distributions()
lab_b.intensity_distributions()
lab_c.intensity_distributions()
#%%
local_dasen_paramA = lab_a.local_normalisation_parameters()
local_dasen_paramB = lab_b.local_normalisation_parameters()
local_dasen_paramC = lab_c.local_normalisation_parameters()
#%% md
# Server side
#%%
probe_type_means = serv.aggregate_QN_means(local_dasen_paramA, local_dasen_paramB, local_dasen_paramC)
#%% md
# Client side
#%%
lab_a.final_normalisation(probe_type_means)
lab_b.final_normalisation(probe_type_means)
lab_c.final_normalisation(probe_type_means)
#%%
# save the betas for testing
lab_a.betas.to_csv(os.path.join(output, f"{split_type}_split1_betas.csv"))
lab_b.betas.to_csv(os.path.join(output, f"{split_type}_split2_betas.csv"))
lab_c.betas.to_csv(os.path.join(output, f"{split_type}_split3_betas.csv"))
#%% md
# ## EWAS - Linear regression model
#%% md
# Client side
#%%
lab_a.local_xtx_xty()
lab_b.local_xtx_xty()
lab_c.local_xtx_xty()
#%% md
# Server side
#%%
serv.global_regression_parameter((lab_a.xtx, lab_a.xty), (lab_b.xtx, lab_b.xty), (lab_c.xtx, lab_c.xty))
#%% md
# #### Client side
#%% md
# calculate the local sse and covariance of the regression coefficients
#%%
SSE_a,cov_coef_a = lab_a.compute_SSE_and_cov_coef(serv.beta)
SSE_b,cov_coef_b = lab_b.compute_SSE_and_cov_coef(serv.beta)
SSE_c,cov_coef_c = lab_c.compute_SSE_and_cov_coef(serv.beta)
SSE_list = [SSE_a, SSE_b, SSE_c]
cov_coef_list = [cov_coef_a, cov_coef_b, cov_coef_c]
#%% md
# #### Server side
#%% md
# calculate the global SSE and covariance of the regression coefficients
#%%
serv.aggregate_SSE_and_cov_coef(SSE_list,cov_coef_list)
#%%
np.savetxt(os.path.join(output, f"{split_type}_splits_model_matrix.csv"), serv.beta, delimiter=",")
np.savetxt(os.path.join(output, f"{split_type}_splits_model_matrix_xty.csv"), serv.global_xty, delimiter=",")
#%% md
# Make and fit the contrasts to the linear model
#%%
contrasts_mat = serv.make_contrasts(contrasts=[(["AD"],["CTRL"])])
serv.fit_contasts(contrasts_mat.values)
#%%
contrasts_mat.to_csv(os.path.join(output, f"{split_type}_splits_contrastmat.csv"))
np.savetxt(os.path.join(output, f"{split_type}_splits_model_matrix_contractfit.csv"), serv.beta, delimiter=",")
np.savetxt(os.path.join(output, f"{split_type}_splits_model_matrix_xty_contractfit.csv"), serv.global_xty, delimiter=",")
#%% md
# Calculate the P-values
#%%
serv.eBayes()
#%% md
# Get the results table
#%%
serv.table 
#%%
serv.table.to_csv(os.path.join(output, f"{split_type}_splits_EWAS_results.csv"))