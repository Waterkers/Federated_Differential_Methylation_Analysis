#%%
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
#%% md
# ## Evaluation of the workflows
# The performance of the python and federated implmentation of the workflow provided by the x group at Exeter university are compared to the original r-version of this workflow.
# Additionally the performance of the federated workflow is tested in several cases of sample size and class label imbalance.
#%%
# define helper functions
def groundTruthNew(R_EWAS_results, diff_meth_thresh=0.05):
    probe_difference = pd.Series.abs(R_EWAS_results.loc[:, "logFC"])
    DM_probes = probe_difference[probe_difference >= diff_meth_thresh].index.values
    DM_probes = set(DM_probes).intersection(
        set(R_EWAS_results.loc[R_EWAS_results["adj.P.Val"] <= 0.05, :].index.values))
    notDM_probes = set(probe_difference[probe_difference >= diff_meth_thresh].index.values).difference(DM_probes)
    return {"DM": DM_probes, "Not DM": notDM_probes}


def EWAS_metrics_p(EWAS_results, ground_truth, diff_meth_thresh=0.05, multiIndex=False):
    # control_probe_means = Normalised_betas.loc[:, phenotype.index[phenotype["Diagnosis"] == " CTRL"]].mean(axis=1)
    # AD_probe_means = Normalised_betas.loc[:, phenotype.index[phenotype["Diagnosis"] != " CTRL"]].mean(axis=1)

    probe_difference = pd.Series.abs(EWAS_results.loc[:, "logFC"])
    EWAS_sig = EWAS_results.loc[EWAS_results["adj.P.Val"] <= 0.05, :].index

    good_probes = probe_difference[probe_difference >= diff_meth_thresh].index
    good_probes = set(good_probes).intersection(set(EWAS_sig))
    not_good = set(probe_difference[probe_difference >= diff_meth_thresh].index.values).difference(good_probes)

    # true positives - EWAS significant and differentially methylated in python/fed and in original R
    true_positive = len(ground_truth["DM"].intersection(good_probes))
    # true negatives - EWAS not significant and not differentially methylated in python/fed and in original R
    true_negative = len(ground_truth["Not DM"].intersection(not_good))
    # false positives - EWAS significant and not differentially methylated in python/fed and in original R
    false_positve = len(ground_truth["Not DM"].intersection(good_probes))
    # false negatives = EWAS not significant and differentially methylated in python/fed and in original R
    false_negative = len(ground_truth["DM"].intersection(not_good))

    if (true_positive + true_negative + false_positve + false_negative) > 0:
        Acc = (true_positive + true_negative) / (true_positive + true_negative + false_positve + false_negative)
    else:
        Acc = 0
    if (true_positive + false_positve) > 0:
        Pre = true_positive / (true_positive + false_positve)
    else:
        Pre = 0
    if (true_positive + false_negative) > 0:
        Rec = true_positive / (true_positive + false_negative)
    else:
        Rec = 0
    if Pre and Rec:
        F1 = 2 * (Rec * Pre) / (Rec + Pre)
    else:
        F1 = 0

    return {"TP": true_positive, "TN": true_negative, "FP": false_positve, "FN": false_negative, "Accuracy": Acc,
            "Precision": Pre, "Recall": Rec, "F1": F1}

def plot_pval_distributions(experiment_data:pd.DataFrame, ground_truth_data:pd.DataFrame, identifier:str, distortion:str,
                          save_file_name:str=None, return_axis:bool=False):
    #scatterplot p-values distributions r and fed python in the different scenarios
    fig, (ax1, ax2) = plt.subplots(ncols = 1, nrows=2, figsize=(5,7.5))
    # ground truth
    ax1.scatter(x=ground_truth_data.index.values, y = ground_truth_data.loc[:, "adj.P.Val"], c="r")
    ax1.set_ylim((0.00,1.00))
    ax1.set_ylabel("Adjusted P-value")
    ax1.set_xlabel("Probes")
    ax1.set_xticks([])
    ax1.set_title(f"{identifier} - Ground truth")
    # distorted - fed
    ax2.scatter(x=experiment_data.index.values, y = experiment_data.loc[:, "adj.P.Val"], c="b")
    ax2.set_ylim((0.00,1.00))
    ax2.set_ylabel("Adjusted P-value")
    ax2.set_xlabel("Probes")
    ax2.set_xticks([])
    ax2.set_title(f"{distortion} - Federated")
    if save_file_name:
        fig.savefig(save_file_name)
        plt.close()
    if return_axis:
        return ax1, ax2

def plot_pval_probes_diff_gt_fed(experiment_data:pd.DataFrame, ground_truth_data:pd.DataFrame,
                            identifier:str, distortion:str, alpha:float=0.05, save_file_name:str=None, return_axis:bool=False):
    diff_ids = set(
        ground_truth_data.loc[ground_truth_data["adj.P.Val"] <= alpha, :].index.values).difference(
        set(experiment_data.loc[experiment_data["adj.P.Val"] <= alpha, :].index.values))
    #scatterplot with p-values for the diff ids in R and fed-python
    fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, figsize=(5,7.5))
    # ground truth
    ax1.scatter(x=list(diff_ids), y = ground_truth_data.loc[diff_ids, "adj.P.Val"], c="r")
    ax1.set_ylim((0.00,0.15))
    ax1.set_ylabel("Adjusted P-value")
    ax1.set_xlabel("Probes")
    ax1.set_xticks([])
    ax1.set_title(f"{identifier} - Ground truth")
    # distortion - fed
    ax2.scatter(x=list(diff_ids), y = experiment_data.loc[diff_ids, "adj.P.Val"], c="b")
    ax2.set_ylim((0.00,0.15))
    ax2.set_ylabel("Adjusted P-value")
    ax2.set_xlabel("Probes")
    ax2.set_xticks([])
    ax2.set_title(f"{distortion} - Federated")
    if save_file_name:
        fig.savefig(save_file_name)
        plt.close()
    if return_axis:
        return ax1, ax2

def calculate_metrics_dataset(experiment_data, ground_truth_data, identifier, save_file_name, sheet_name):
    thresholds = {"0.01": 0.01, "0.05": 0.05, "0.1": 0.1, "0.15": 0.15, "0.2": 0.20}
    inner = {}
    for i in thresholds:
        ground_truth = groundTruthNew(ground_truth_data, diff_meth_thresh=thresholds[i])
        metrics = EWAS_metrics_p(experiment_data, diff_meth_thresh=thresholds[i], ground_truth=ground_truth,
                                 multiIndex=True)
        inner[i] = metrics
    result = pd.DataFrame.from_dict(inner).T
    # save the metric tables to excel
    with pd.ExcelWriter(save_file_name,
                        mode="w", engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name=(identifier + sheet_name))
    return result



def create_performance_metrics_excel(experimental_data_dict:dict, ground_truth_data:pd.DataFrame, save_file_name:str, return_results:bool=True):
    results = {}
    for data in experimental_data_dict:
        results[data] = calculate_metrics_dataset(datasets[data], ground_truth_data, data,
                                                  save_file_name,
                                                  "EWASMetrics")

    # combine the metrics results into one dataframe
    performance_results = pd.concat(list(results.values()), keys=list(results.keys()))
    with pd.ExcelWriter(save_file_name, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer:
        performance_results.to_excel(writer, sheet_name='combinedPerformanceResults')
    if return_results:
        return performance_results
#%%
#define helper variables
phenotype_recoding = {'GSE66351':{'Diagnosis':{"diagnosis: AD":"AD",
                                               "diagnosis: CTRL":"CTRL"},
                                  'Sex':{"":"F",
                                         "":"M"}},
                      "GSE105109":{'Diagnosis':{"post-mortem diagnosis: Alzheimer's disease":"AD",
                                                "post-mortem diagnosis: Control":"CTRL"}},
                      'GSE134379':{"Diagnosis":{"diagnosis: AD":"AD",
                                                "diagnosis: ND":"CTRL"}}}


#%%
# compile all the steps that need to be performed per dataset on ONE first
# then create function around it to make it easy to call on all of them
# GSE66351_pheno = pd.read_csv("E:\\Msc Systems Biology\\MSB5000_Master_Thesis\\Practical work\\Data\\Data_Full_Datasets\\GSE66351\\Reduced_Pheno_Info.csv", index_col=0)
# recoded_pheno = GSE66351_pheno.copy(deep=True)
# for column in phenotype_recoding['GSE66351']:
#     recoded_pheno[column] = recoded_pheno[column].map(phenotype_recoding['GSE66351'][column])

# make sure all the phenotype files used the same codes for the diagnosis

#%%
# read in the results from the different implementations
original_r_GSE66351_E = pd.read_csv("/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/Small_Results_dataset.csv", index_col=0)
# python
python_central_GSE66351_E = pd.read_csv("/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_eBayesTopTableResult.csv", index_col=0)
distortions = ['balanced', 'strong', 'mild']
federated_balanced = pd.read_csv("/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/balanced_splits_EWAS_results.csv", index_col=0)
#federated_mild_imbalance = pd.read_csv("E:\Msc Systems Biology\MSB5000_Master_Thesis\Practical work\Data\GSE66351_Fed\mild_splits_EWAS_results.csv", index_col=0)
#federated_strong_imbalance = pd.read_csv("E:\Msc Systems Biology\MSB5000_Master_Thesis\Practical work\Data\GSE66351_Fed\strong_splits_EWAS_results.csv", index_col=0)
#%%
# calculate diff and non-dfff methylated probes for all of the methods for checking purposes
GSE66351_groundTruth = groundTruthNew(original_r_GSE66351_E)
python_ground_truth = groundTruthNew(python_central_GSE66351_E)
federatedGroundTruth = groundTruthNew(federated_balanced)
#%%
# loop through different methylation cut-offs to compare perfomance of R and fed
datasets = {"Python":python_central_GSE66351_E,"Balanced":federated_balanced} #, "Mild Imbalance":federated_mild_imbalance, "Strong Imbalance":federated_strong_imbalance}
rGroundTruth_results = create_performance_metrics_excel(datasets, original_r_GSE66351_E,
                                 "/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_rGroundTruth_EvaluationResults.xlsx",
                                 return_results=True)
# compare against R central run
for data in datasets:
    plot_pval_distributions(datasets[data], original_r_GSE66351_E, 'GSE66351', data,
                            "/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_rGroundTruth_pvalDistributions.png")
    plot_pval_probes_diff_gt_fed(datasets[data], original_r_GSE66351_E, 'GSE66351', data,
                                 save_file_name="/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_rGroundTruth_pvalProbeDifferences.png")
#%%
datasets = {"R":original_r_GSE66351_E,"Balanced":federated_balanced}
pythonGroundTruth_results = create_performance_metrics_excel(datasets, python_central_GSE66351_E,
                                 "/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_pythonGroundTruth_EvaluationResults.xlsx",
                                 return_results=True)
# compare against python central run
for data in datasets:
    plot_pval_distributions(datasets[data], python_central_GSE66351_E, 'GSE66351', data,
                            "/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_pythonGroundTruth_pvalDistributions.png")
    plot_pval_probes_diff_gt_fed(datasets[data], python_central_GSE66351_E, 'GSE66351', data,
                                 save_file_name="/home/silke/Documents/Fed_EWAS/Data/EWASResultsServer/GSE66351_noCohortEffects_pythonGroundTruth_pvalProbeDifferences.png")
#%%