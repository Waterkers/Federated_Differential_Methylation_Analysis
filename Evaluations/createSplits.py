import pandas as pd
import os
from random import sample, seed

def generateEvenSplits(phenoPath:str,
                       methFilteredPath:str,
                       unMethFilteredPath:str,
                       betaFilteredPath:str,
                       designPath:str,
                       designLocalPath:str,
                       fullDesignPath:str,
                       fullDesignLocalPath:str,
                       outputPath:str,
                       identifier,
                       n_splits:int=3):
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    if n_splits!=3:
        print('For now only 3 splits are supported, setting n_splits=3')
        n_splits=3

    GSE66351_pheno = pd.read_csv(
        phenoPath,
        index_col="Sample_ID")
    GSE66351_splits_pheno = GSE66351_pheno.copy()
    GSE66351_splits_pheno["split"] = "Split_3"

    N_total = GSE66351_splits_pheno.shape[0] * 1.0 * 1.0
    N_ad = GSE66351_splits_pheno.loc[GSE66351_splits_pheno["Diagnosis"] == " AD", :].shape[0] * 1.0 * 1.0
    random_state = 42
    seed(random_state)
    sizes = [1, 1, 1]
    ad_freqs = [0.53, 0.53, 0.53]

    Sizes = []
    n_ad = []
    for i in range(0, n_splits - 1):
        s = int(N_total * sizes[i] / sum(sizes))
        Sizes.append(s)
        n_ad.append(int(s * ad_freqs[i]))

    Sizes.append(int(N_total - sum(Sizes)))
    n_ad.append(int(N_ad - sum(n_ad)))
    print(Sizes, sum(Sizes))
    print(n_ad, sum(n_ad))

    splits = {}
    ad = set(GSE66351_pheno.loc[GSE66351_pheno["Diagnosis"] == "diagnosis: AD", :].index.values)
    other = set(GSE66351_pheno.index.values).difference(ad)  # .difference(fem)
    for i in range(0, n_splits - 1):
        b = set(sample(ad, n_ad[i]))
        ad = ad.difference(b)
        o = set(sample(other, Sizes[i] - n_ad[i]))
        other = other.difference(o)
        sele_samples = b | o
        GSE66351_splits_pheno.loc[sele_samples, "split"] = "Split_" + str(i + 1)
        GSE66351_splits_pheno["Split_" + str(i + 1)] = 0
        GSE66351_splits_pheno.loc[sele_samples, "Split_" + str(i + 1)] = 1
    print(GSE66351_splits_pheno[["split", "Diagnosis"]].groupby("split")["Diagnosis"].value_counts())

    # save the dataset splits
    meth = pd.read_csv(
        methFilteredPath,
        index_col=0)
    unmeth = pd.read_csv(
        unMethFilteredPath,
        index_col=0)
    beta = pd.read_csv(
        betaFilteredPath,
        index_col=0)
    design = pd.read_csv(
        designPath,
        index_col=0)
    design_local = pd.read_csv(
        designLocalPath,
        index_col=0)
    design_full = pd.read_csv(
        fullDesignPath,
        index_col=0)
    design_full_local = pd.read_csv(
        fullDesignLocalPath,
        index_col=0)

    output_dir = os.path.join(outputPath, f"{identifier}_splits")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for i in range(n_splits):
        s = "Split_" + str(i + 1)
        GSE66351_pheno = GSE66351_splits_pheno.loc[GSE66351_splits_pheno["split"] == s, :]
        samples = sorted(GSE66351_pheno.index.values)
        print("these samples: %s are included in split %s" % (len(samples), i + 1))
        GSE66351_splits_pheno.loc[samples, :].to_csv(output_dir + "/" + s + "_pheno.csv")
        meth.loc[:, samples].to_csv(output_dir + "/" + s + "_methylated.csv")
        unmeth.loc[:, samples].to_csv(output_dir + "/" + s + "_unmethylated.csv")
        beta.loc[:, samples].to_csv(output_dir + "/" + s + "_betas.csv")
        design.loc[samples, :].to_csv(output_dir + "/" + s + "_design.csv")
        design_full.loc[samples, :].to_csv(output_dir + "/" + s + "_Full_design.csv")

    central_design = design_local.copy()
    central_design["Cohort_effect"] = GSE66351_splits_pheno["split"]
    central_design["Split1"] = 0
    central_design.loc[central_design["Cohort_effect"] == "Split_1", "Split1"] = 1
    central_design["Split2"] = 0
    central_design.loc[central_design["Cohort_effect"] == "Split_2", "Split2"] = 1
    central_design.drop(columns="Cohort_effect", inplace=True)
    central_design.to_csv(os.path.join(output_dir + "/" + "central_design_matrix.csv"))

    central_full_design = design_full_local.copy()
    central_full_design["Cohort_effect"] = GSE66351_splits_pheno["split"]
    central_full_design["Split1"] = 0
    central_full_design.loc[central_full_design["Cohort_effect"] == "Split_1", "Split1"] = 1
    central_full_design["Split2"] = 0
    central_full_design.loc[central_full_design["Cohort_effect"] == "Split_2", "Split2"] = 1
    central_full_design.drop(columns="Cohort_effect", inplace=True)
    central_full_design.to_csv(os.path.join(output_dir + "/" + "full_central_design_matrix.csv"))