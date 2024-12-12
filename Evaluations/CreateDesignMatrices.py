import pandas as pd
import os


def createDesignMatrix66351(pheno_df_path: str, small: bool = True, federated: bool = False, per_region: bool = False,
                       output_path: str = None):
    pheno = pd.read_csv(pheno_df_path, index_col="Sample_ID", low_memory=False)
    pheno["Diagnosis"] = pheno.loc[:, "Diagnosis"].str.strip()
    pheno["Sex"] = pheno.loc[:, "Sex"].str.strip()
    pheno["Brain_region"] = pheno.loc[:, "Brain_region"].str.strip()
    pheno["Brain_region"] = pheno.loc[:, "Brain_region"].str.replace(" ", "")
    x = pheno.loc[:, ["Diagnosis", "Age", "Sex", "sentrix_id",
                      "Brain_region"]]  # design matrix with the dependent/explainatory variables to be included in the model
    # The design matrix needs to consist of numeric representations of the covariates to be included in the model, i.e. binary diagnosis, binary sex, dummy sentrix etc.
    x["AD"] = 0
    x.loc[x["Diagnosis"] == "AD", "AD"] = 1  # create binary diagnosis with 1 = AD and 0 = CTR
    x["CTRL"] = 0
    x.loc[x["Diagnosis"] == "CTRL", "CTRL"] = 1
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

if __name__ == "__main__":
    createDesignMatrix66351(pheno_df_path='/home/silke/Documents/Fed_EWAS/Data/QC_GSE66351_half/Reduced_Pheno_Info.csv',
                        small=True, federated=False, per_region=False,
                        output_path='/home/silke/Documents/Fed_EWAS/Data/QC_GSE66351_half')