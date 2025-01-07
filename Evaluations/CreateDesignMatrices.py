import pandas as pd
import os
import argparse

def createDesignMatrix66351(pheno_df_path: str, small: bool = True, federated: bool = False, per_region: bool = False,
                       output_path: str = None, half: bool = False):
    #if half:
    pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False)
    #else:
        #pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False, sep='\t').T
    pheno["Diagnosis"] = pheno.loc[:, "Diagnosis"].str.strip()
    pheno["Sex"] = pheno.loc[:, "Sex"].str.strip()
    pheno["Brain_region"] = pheno.loc[:, "Brain_region"].str.strip()
    pheno["Brain_region"] = pheno.loc[:, "Brain_region"].str.replace(" ", "")
    x = pheno.loc[:, ["Diagnosis", "Age", "Sex", "sentrix_id",
                      "Brain_region"]].copy(deep=True)  # design matrix with the dependent/explainatory variables to be included in the model
    # The design matrix needs to consist of numeric representations of the covariates to be included in the model, i.e. binary diagnosis, binary sex, dummy sentrix etc.
    x["AD"] = 0
    x.loc[x["Diagnosis"] == "diagnosis: AD", "AD"] = 1  # create binary diagnosis with 1 = AD and 0 = CTR
    x["CTRL"] = 0
    x.loc[x["Diagnosis"] == "diagnosis: CTRL", "CTRL"] = 1
    x["Sex"] = x.loc[:, "Sex"].replace("^[^:]*:", "", regex=True)
    x.loc[x["Sex"] == " F", "Sex"] = 1
    x.loc[x["Sex"] == " M", "Sex"] = 0
    # x.loc[x["Sex"] == " F", "Sex"] = 1 #create binary sex with 1 = F and 0 = M

    # turn the age variable into a continuous numerical variable without any leftover text
    x["Age"] = x.loc[:, "Age"].replace("^[^:]*:", "", regex=True)
    x["Age"] = pd.to_numeric(x["Age"])
    tissues = pheno["Brain_region"].unique()
    if small:
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Brain_region"] == (" " + tis), :]
                    t = x.drop(columns=["Brain_region"])
                    if output_path:
                        # use output path
                        t.to_csv(os.path.join(output_path, "Small_%s_EWAS_design.csv" % (tis)))
                    else:
                        # otherwise save in current working dir
                        t.to_csv("Small_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    # use output path
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "sentrix_id"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "sentrix_id"]].to_csv("Small_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_Tcen = x.copy()
                    x_Tcen = x_Tcen.loc[pheno["Brain_region"] == tis, :]
                    unique_ids = x_Tcen["sentrix_id"].unique()
                    for id in unique_ids[:-1]:
                        x_Tcen[id] = (x_Tcen["sentrix_id"] == id).astype(int)
                    t_cen = x_Tcen.drop(columns=["Diagnosis", "Brain_region", "sentrix_id"])
                    if output_path:
                        t_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_cen.to_csv("Small_%s_EWAS_design_local.csv" % (tis))
            else:
                x_cen = x.copy()
                unique_ids = x_cen["sentrix_id"].unique()
                include = unique_ids[:-1]
                for id in unique_ids[:-1]:
                    x_cen[id] = (x_cen["sentrix_id"] == id).astype(int)
                x_cen.drop(columns=["Diagnosis", "sentrix_id", "Brain_region"], inplace=True)
                if output_path:
                    x_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_cen.to_csv("Small_EWAS_design_local.csv")
    else:
        # now on to the full design matrix
        x_large = x.copy()
        x_large["Cell_Type1"] = pd.to_numeric(pheno["Cell_Type.CT1"])
        x_large["Cell_Type2"] = pd.to_numeric(pheno["Cell_Type.CT2"])
        x_large["Cell_Type3"] = pd.to_numeric(pheno["Cell_Type.CT3"])
        x_large["Cell_Type4"] = pd.to_numeric(pheno["Cell_Type.CT4"])
        x_large["Cell_Type5"] = pd.to_numeric(pheno["Cell_Type.CT5"])
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Brain_region"] == (" " + tis), :]
                    t_large = x_large.drop(columns=["Brain_region"])
                    if output_path:
                        t_large.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_large.to_csv("Full_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "sentrix_id",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "sentrix_id",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        "Full_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_large_Tcen = x_large.copy()
                    x_large_Tcen = x_large_Tcen.loc[pheno["Brain_region"] == tis, :]
                    unique_ids = x_large_Tcen["sentrix_id"].unique()
                    for id in unique_ids[:-1]:
                        x_large_Tcen[id] = (x_large_Tcen["sentrix_id"] == id).astype(int)
                    t_large_cen = x_large_Tcen.drop(columns=["Diagnosis", "Brain_region", "sentrix_id"])
                    if output_path:
                        t_large_cen.to_csv(os.path.join(output_path, "Full_%s_EWAS_design_local.csv" % (tis)))
                    else:
                        t_large_cen.to_csv("Full_%s_EWAS_design_local.csv" % (tis))
            else:
                x_large_cen = x_large.copy()
                unique_ids = x_large_cen["sentrix_id"].unique()
                for id in unique_ids[:-1]:
                    x_large_cen[id] = (x_large_cen["sentrix_id"] == id).astype(int)
                x_large_cen.drop(columns=["Diagnosis", "sentrix_id", "Brain_region"], inplace=True)
                if output_path:
                    x_large_cen.to_csv(os.path.join(output_path, "Full_EWAS_design_local.csv"))
                else:
                    x_large_cen.to_csv("Full_EWAS_design_local.csv")

def createDesignMatrix105109(pheno_df_path: str, small: bool = True, federated: bool = False, per_region: bool = False,
                       output_path: str = None, half: bool = False):
    # if half:
    pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False)
    # else:
    # pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False, sep='\t').T
    pheno["Diagnosis"] = pheno.loc[:, "Diagnosis"].str.strip()
    pheno["Sex"] = pheno.loc[:, "Sex"].str.strip()
    pheno["Source_Tissue"] = pheno.loc[:, "Source_Tissue"].str.strip()
    pheno["Source_Tissue"] = pheno.loc[:, "Source_Tissue"].str.replace(" ", "")
    x = pheno.loc[:, ["Diagnosis", "Age", "Sex", "Sentrix_ID",
                      "Source_Tissue"]]  # design matrix with the dependent/explainatory variables to be included in the model
    # The design matrix needs to consist of numeric representations of the covariates to be included in the model, i.e. binary diagnosis, binary sex, dummy sentrix etc.
    x["AD"] = 0
    x.loc[x["Diagnosis"] == "post-mortem diagnosis: Alzheimer's disease", "AD"] = 1  # create binary diagnosis with 1 = AD and 0 = CTR
    x["CTRL"] = 0
    x.loc[x["Diagnosis"] == "post-mortem diagnosis: Control", "CTRL"] = 1
    x.loc[x["Sex"] == "gender: F", "Sex"] = 1
    x.loc[x["Sex"] == "gender: M", "Sex"] = 0
    # x.loc[x["Sex"] == " F", "Sex"] = 1 #create binary sex with 1 = F and 0 = M

    # turn the age variable into a continuous numerical variable without any leftover text
    x["Age"] = x.loc[:, "Age"].replace("^[^:]*:", "", regex=True)
    x["Age"] = pd.to_numeric(x["Age"])
    tissues = pheno["Source_Tissue"].unique()
    if small:
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Source_Tissue"] == (" " + tis), :]
                    t = x.drop(columns=["Source_Tissue"])
                    if output_path:
                        # use output path
                        t.to_csv(os.path.join(output_path, "Small_%s_EWAS_design.csv" % (tis)))
                    else:
                        # otherwise save in current working dir
                        t.to_csv("Small_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    # use output path
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID"]].to_csv("Small_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_Tcen = x.copy()
                    x_Tcen = x_Tcen.loc[pheno["Source_Tissue"] == tis, :]
                    unique_ids = x_Tcen["Sentrix_ID"].unique()
                    for id in unique_ids[:-1]:
                        x_Tcen[id] = (x_Tcen["Sentrix_ID"] == id).astype(int)
                    t_cen = x_Tcen.drop(columns=["Diagnosis", "Source_Tissue", "Sentrix_ID"])
                    if output_path:
                        t_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_cen.to_csv("Small_%s_EWAS_design_local.csv" % (tis))
            else:
                x_cen = x.copy()
                unique_ids = x_cen["Sentrix_ID"].unique()
                include = unique_ids[:-1]
                for id in unique_ids[:-1]:
                    x_cen[id] = (x_cen["Sentrix_ID"] == id).astype(int)
                x_cen.drop(columns=["Diagnosis", "Sentrix_ID", "Source_Tissue"], inplace=True)
                if output_path:
                    x_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_cen.to_csv("Small_EWAS_design_local.csv")
    else:
        # now on to the full design matrix
        x_large = x.copy()
        x_large["Cell_Type1"] = pd.to_numeric(pheno["Cell_Type.CT1"])
        x_large["Cell_Type2"] = pd.to_numeric(pheno["Cell_Type.CT2"])
        x_large["Cell_Type3"] = pd.to_numeric(pheno["Cell_Type.CT3"])
        x_large["Cell_Type4"] = pd.to_numeric(pheno["Cell_Type.CT4"])
        x_large["Cell_Type5"] = pd.to_numeric(pheno["Cell_Type.CT5"])
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Source_Tissue"] == (" " + tis), :]
                    t_large = x_large.drop(columns=["Source_Tissue"])
                    if output_path:
                        t_large.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_large.to_csv("Full_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        "Full_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_large_Tcen = x_large.copy()
                    x_large_Tcen = x_large_Tcen.loc[pheno["Source_Tissue"] == tis, :]
                    unique_ids = x_large_Tcen["Sentrix_ID"].unique()
                    for id in unique_ids[:-1]:
                        x_large_Tcen[id] = (x_large_Tcen["Sentrix_ID"] == id).astype(int)
                    t_large_cen = x_large_Tcen.drop(columns=["Diagnosis", "Source_Tissue", "Sentrix_ID"])
                    if output_path:
                        t_large_cen.to_csv(os.path.join(output_path, "Full_%s_EWAS_design_local.csv" % (tis)))
                    else:
                        t_large_cen.to_csv("Full_%s_EWAS_design_local.csv" % (tis))
            else:
                x_large_cen = x_large.copy()
                unique_ids = x_large_cen["Sentrix_ID"].unique()
                for id in unique_ids[:-1]:
                    x_large_cen[id] = (x_large_cen["Sentrix_ID"] == id).astype(int)
                x_large_cen.drop(columns=["Diagnosis", "Sentrix_ID", "Source_Tissue"], inplace=True)
                if output_path:
                    x_large_cen.to_csv(os.path.join(output_path, "Full_EWAS_design_local.csv"))
                else:
                    x_large_cen.to_csv("Full_EWAS_design_local.csv")

def createDesignMatrix134379(pheno_df_path: str, small: bool = True, federated: bool = False, per_region: bool = False,
                       output_path: str = None, half: bool = False):
    # if half:
    pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False)
    # else:
    # pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False, sep='\t').T
    pheno["Diagnosis"] = pheno.loc[:, "Diagnosis"].str.strip()
    pheno["Sex"] = pheno.loc[:, "Sex"].str.strip()
    pheno["Source_region"] = pheno.loc[:, "Source_region"].str.strip()
    pheno["Source_region"] = pheno.loc[:, "Source_region"].str.replace(" ", "")
    x = pheno.loc[:, ["Diagnosis", "Age", "Sex", "Sentrix_ID",
                      "Source_region"]]  # design matrix with the dependent/explainatory variables to be included in the model
    # The design matrix needs to consist of numeric representations of the covariates to be included in the model, i.e. binary diagnosis, binary sex, dummy sentrix etc.
    x["AD"] = 0
    x.loc[x["Diagnosis"] == "diagnosis: AD", "AD"] = 1  # create binary diagnosis with 1 = AD and 0 = CTR
    x["CTRL"] = 0
    x.loc[x["Diagnosis"] == "diagnosis: ND", "CTRL"] = 1
    x.loc[x["Sex"] == "gender: F", "Sex"] = 1
    x.loc[x["Sex"] == "gender: M", "Sex"] = 0
    # x.loc[x["Sex"] == " F", "Sex"] = 1 #create binary sex with 1 = F and 0 = M

    # turn the age variable into a continuous numerical variable without any leftover text
    x["Age"] = x.loc[:, "Age"].replace("^[^:]*:", "", regex=True)
    x["Age"] = pd.to_numeric(x["Age"])
    tissues = pheno["Source_region"].unique()
    if small:
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Source_region"] == (" " + tis), :]
                    t = x.drop(columns=["Source_region"])
                    if output_path:
                        # use output path
                        t.to_csv(os.path.join(output_path, "Small_%s_EWAS_design.csv" % (tis)))
                    else:
                        # otherwise save in current working dir
                        t.to_csv("Small_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    # use output path
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID"]].to_csv("Small_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_Tcen = x.copy()
                    x_Tcen = x_Tcen.loc[pheno["Source_region"] == tis, :]
                    unique_ids = x_Tcen["Sentrix_ID"].unique()
                    for id in unique_ids[:-1]:
                        x_Tcen[id] = (x_Tcen["Sentrix_ID"] == id).astype(int)
                    t_cen = x_Tcen.drop(columns=["Diagnosis", "Source_region", "Sentrix_ID"])
                    if output_path:
                        t_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_cen.to_csv("Small_%s_EWAS_design_local.csv" % (tis))
            else:
                x_cen = x.copy()
                unique_ids = x_cen["Sentrix_ID"].unique()
                include = unique_ids[:-1]
                for id in unique_ids[:-1]:
                    x_cen[id] = (x_cen["Sentrix_ID"] == id).astype(int)
                x_cen.drop(columns=["Diagnosis", "Sentrix_ID", "Source_region"], inplace=True)
                if output_path:
                    x_cen.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_cen.to_csv("Small_EWAS_design_local.csv")
    else:
        # now on to the full design matrix
        x_large = x.copy()
        x_large["Cell_Type1"] = pd.to_numeric(pheno["Cell_Type.CT1"])
        x_large["Cell_Type2"] = pd.to_numeric(pheno["Cell_Type.CT2"])
        x_large["Cell_Type3"] = pd.to_numeric(pheno["Cell_Type.CT3"])
        x_large["Cell_Type4"] = pd.to_numeric(pheno["Cell_Type.CT4"])
        x_large["Cell_Type5"] = pd.to_numeric(pheno["Cell_Type.CT5"])
        if federated:
            if per_region:
                for tis in tissues:
                    # x.loc[pheno["Source_region"] == (" " + tis), :]
                    t_large = x_large.drop(columns=["Source_region"])
                    if output_path:
                        t_large.to_csv(os.path.join(output_path, "Small_EWAS_design.csv"))
                    else:
                        t_large.to_csv("Full_%s_EWAS_design.csv" % (tis))

            else:
                if output_path:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        os.path.join(output_path, "Small_EWAS_design.csv"))
                else:
                    x_large.loc[:, ["AD", "CTRL", "Age", "Sex", "Sentrix_ID",
                                    "Cell_Type1", "Cell_Type2", "Cell_Type3", "Cell_Type4", "Cell_Type5"]].to_csv(
                        "Full_EWAS_design.csv")
        else:
            if per_region:
                for tis in tissues:
                    x_large_Tcen = x_large.copy()
                    x_large_Tcen = x_large_Tcen.loc[pheno["Source_region"] == tis, :]
                    unique_ids = x_large_Tcen["Sentrix_ID"].unique()
                    for id in unique_ids[:-1]:
                        x_large_Tcen[id] = (x_large_Tcen["Sentrix_ID"] == id).astype(int)
                    t_large_cen = x_large_Tcen.drop(columns=["Diagnosis", "Source_region", "Sentrix_ID"])
                    if output_path:
                        t_large_cen.to_csv(os.path.join(output_path, "Full_%s_EWAS_design_local.csv" % (tis)))
                    else:
                        t_large_cen.to_csv("Full_%s_EWAS_design_local.csv" % (tis))
            else:
                x_large_cen = x_large.copy()
                unique_ids = x_large_cen["Sentrix_ID"].unique()
                for id in unique_ids[:-1]:
                    x_large_cen[id] = (x_large_cen["Sentrix_ID"] == id).astype(int)
                x_large_cen.drop(columns=["Diagnosis", "Sentrix_ID", "Source_region"], inplace=True)
                if output_path:
                    x_large_cen.to_csv(os.path.join(output_path, "Full_EWAS_design_local.csv"))
                else:
                    x_large_cen.to_csv("Full_EWAS_design_local.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pheno_df_path', type=str,help='path to the Pheno_Info.csv file')
    parser.add_argument('output_path', type=str,help='path to folder to store the design matrices')
    parser.add_argument('-s', '--small', action='store_true', help='create small design matrices')
    parser.add_argument('-f', '--federated', action='store_true', help='create federated design matrices')
    parser.add_argument('--per_region', action='store_true', help='create per region design matrices')
    args = parser.parse_args()
    if '66351' in args.pheno_df_path:
        createDesignMatrix66351(pheno_df_path=args.pheno_df_path,
                            small=args.small, federated=args.federated, per_region=args.per_region,
                            output_path=args.output_path)
    elif '105109' in args.pheno_df_path:
        createDesignMatrix105109(pheno_df_path=args.pheno_df_path,
                                small=args.small, federated=args.federated, per_region=args.per_region,
                                output_path=args.output_path)
    elif '134379' in args.pheno_df_path:
        createDesignMatrix134379(pheno_df_path=args.pheno_df_path,
                                small=args.small, federated=args.federated, per_region=args.per_region,
                                output_path=args.output_path)
    else:
        print(f'There is no function to create design matrices for dataset {args.pheno_df_path.split("/")[-1].split("_")[1]}')