#%%
import numpy as np
import os

import pandas as pd

# try:
#     from Federated_Differential_Methylation_Analysis.Central_preprocessing.python_preprocessing_splits import output_dir
# except ModuleNotFoundError:
#     sys.path.append("/cosybio/project/vanElferen/FedEWAS")
#     from Federated_Differential_Methylation_Analysis.Central_preprocessing.python_preprocessing_splits import output_dir
from server import Server
import argparse
from client import Client
#%%
#import pyximport
#pyximport.install(setup_args={"script_args" : ["--verbose"]})
#from linbinR import fast_linbin
#%%

#TODO wrap the whole thing into a function
#TODO include logging or progress monitoring of some kind
#TODO remove all return statements from the client/server functions unless the returned objects are saved to a file for logging
def run_federated_prototype(split_directory:str,
                            output_directory:str,
                            probe_annotation_path:str,
                            cohort_prefix:str='Lab',
                            split:str='balanced',
                            pre_normalised:bool=False,
                            cohort_effected:bool=False,):
    #%% md
    # ## Initialising the clients
    #%%
    # create client
    lab_a = Client(f"{cohort_prefix}_A", os.path.join(split_directory, "Split_1_design.csv"), os.path.join(split_directory, "Split_1_methylated.csv"), os.path.join(split_directory, "Split_1_unmethylated.csv"), probe_annotation_path)
    lab_b = Client(f"{cohort_prefix}_B", os.path.join(split_directory, "Split_2_design.csv"), os.path.join(split_directory, "Split_2_methylated.csv"), os.path.join(split_directory, "Split_2_unmethylated.csv"), probe_annotation_path)
    lab_c = Client(f"{cohort_prefix}_C", os.path.join(split_directory, "Split_3_design.csv"), os.path.join(split_directory, "Split_3_methylated.csv"), os.path.join(split_directory, "Split_3_unmethylated.csv"), probe_annotation_path)
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
    if cohort_effected:
        cohort_effect = serv.find_cohort_effects()
    #%%
    #check client input
    lab_a.input_validation(global_conditions, global_probes)
    lab_b.input_validation(global_conditions, global_probes)
    lab_c.input_validation(global_conditions, global_probes)
    #%%
    if cohort_effected:
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


    if not pre_normalised:
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
        lab_a.betas.to_csv(os.path.join(output, f"{split}_split1{'_central_' if cohort else '_'}betas.csv"))
        lab_b.betas.to_csv(os.path.join(output, f"{split}_split2{'_central_' if cohort else '_'}betas.csv"))
        lab_c.betas.to_csv(os.path.join(output, f"{split}_split3{'_central_' if cohort else '_'}betas.csv"))
    else:
        lab_a.betas = pd.read_csv(os.path.join(split_directory, "Split_1_betas.csv"), index_col=0)
        lab_b.betas = pd.read_csv(os.path.join(split_directory, "Split_2_betas.csv"), index_col=0)
        lab_c.betas = pd.read_csv(os.path.join(split_directory, "Split_3_betas.csv"), index_col=0)
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
    # print statements to figure out why the dimensions between the three splits differ
    print(f'Split 1: xtx {lab_a.xtx.shape}, xty {lab_a.xty.shape}/n '
          f'Split2: xtx {lab_b.xtx.shape}, xty {lab_b.xty.shape}/n '
          f'Split3: xtx {lab_c.xtx.shape}, xty {lab_c.xty.shape}/n ')
    serv.global_regression_parameter((lab_a.xtx, lab_a.xty), (lab_b.xtx, lab_b.xty), (lab_c.xtx, lab_c.xty), cohort=cohort)
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
    np.savetxt(os.path.join(output_directory, f"{split}_splits{'_central_' if cohort else '_'}model_matrix.csv"), serv.beta, delimiter=",")
    np.savetxt(os.path.join(output_directory, f"{split}_splits{'_central_' if cohort else '_'}model_matrix_xty.csv"), serv.global_xty, delimiter=",")
    #%% md
    # Make and fit the contrasts to the linear model
    #%%
    contrasts_mat = serv.make_contrasts(contrasts=[(["AD"],["CTRL"])])
    serv.fit_contasts(contrasts_mat.values)
    #%%
    contrasts_mat.to_csv(os.path.join(output_directory, f"{split}_splits{'_central_' if cohort else '_'}contrastmat.csv"))
    np.savetxt(os.path.join(output_directory, f"{split}_splits{'_central_' if cohort else '_'}model_matrix_contractfit.csv"), serv.beta, delimiter=",")
    np.savetxt(os.path.join(output_directory, f"{split}_splits{'_central_' if cohort else '_'}model_matrix_xty_contractfit.csv"), serv.global_xty, delimiter=",")
    #%% md
    # Calculate the P-values
    #%%
    serv.eBayes()
    #%% md
    # Get the results table
    #%%
    serv.table
    #%%
    serv.table.to_csv(os.path.join(output_directory, f"{split_type}_splits{'_central_' if cohort else '_'}EWAS_results.csv"))

if __name__ == '__main__':
    # args = None
    parser = argparse.ArgumentParser(description='Federated Differential Methylation Analysis')
    parser.add_argument('split_directory', type=str, )
    parser.add_argument('output_dir', type=str, )
    parser.add_argument('probe_annotation_path', type=str, )
    parser.add_argument('split_type', type=str, )
    parser.add_argument('-c', '--cohort', action='store_true')
    parser.add_argument('-n', '--normalised', action='store_true')
    args = parser.parse_args()
    # %%
    if not args:
        split_dir = "/home/silke/Documents/Fed_EWAS/Data/QC_GSE134379_half/GSE134379_splits"
        output = "/home/silke/Documents/Fed_EWAS/Data/QC_GSE134379_half/GSE134379_Fed"
        probeAnnotationPath = "/home/silke/Documents/Fed_EWAS/Data/GSE134379_RAW/GPL13534_HumanMethylation450_15017482_v.1.1.csv"
        split_type = 'balanced'
        cohort = False
        normalised = False
    else:
        split_dir = args.split_directory
        output = args.output_dir
        probeAnnotationPath = args.probe_annotation_path
        split_type = args.split_type
        cohort = args.cohort
        normalised = args.normalised
    # %%
    # check if output directory exists, if not make it
    if not os.path.isdir(output):
        os.mkdir(output)
    run_federated_prototype(split_directory=split_dir,
                            output_directory=output,
                            probe_annotation_path=probeAnnotationPath,
                            split=split_type,
                            cohort_effected=cohort,
                            pre_normalised=normalised,)