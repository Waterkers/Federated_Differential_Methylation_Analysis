import pandas as pd
from random import seed, sample
import os
import argparse
from typing import Union

distortion_ratios = {'GSE66351':{'balanced':{'ad_freqs':[0.53,0.53,0.53],
                                             'sizes':[1,1,1]},
                                 'strong':{'ad_freqs':[0.2,0.3,0.5],
                                             'sizes':[0.15,0.35,0.5]},
                                 'mild':{'ad_freqs':[0.4,0.35,0.25],
                                             'sizes':[0.25,0.35,0.40]}},
                     'GSE134379': {'balanced': {'ad_freqs': [0.55,0.55,0.55],
                                               'sizes': [1, 1, 1]},
                                  'strong': {'ad_freqs': None,
                                             'sizes': None},
                                  'mild': {'ad_freqs': None,
                                           'sizes': None}},
                     'GSE105109': {'balanced': {'ad_freqs': [0.70,0.70,0.70],
                                                'sizes': [1, 1, 1]},
                                   'strong': {'ad_freqs': None,
                                              'sizes': None},
                                   'mild': {'ad_freqs': None,
                                            'sizes': None}},
                     'GSE66351_half': {'balanced': {'ad_freqs': [0.53, 0.53, 0.53],
                                               'sizes': [1, 1, 1]},
                                  'strong': {'ad_freqs': [0.2, 0.3, 0.5],
                                             'sizes': [0.15, 0.35, 0.5]},
                                  'mild': {'ad_freqs': [0.4, 0.35, 0.25],
                                           'sizes': [0.25, 0.35, 0.40]}},
                     'GSE134379_half': {'balanced': {'ad_freqs': [0.55, 0.55, 0.55],
                                                'sizes': [1, 1, 1]},
                                   'strong': {'ad_freqs': None,
                                              'sizes': None},
                                   'mild': {'ad_freqs': None,
                                            'sizes': None}},
                     'GSE105109_half': {'balanced': {'ad_freqs': [0.70, 0.70, 0.70],
                                                'sizes': [1, 1, 1]},
                                   'strong': {'ad_freqs': None,
                                              'sizes': None},
                                   'mild': {'ad_freqs': None,
                                            'sizes': None}},
                     }

label_information = {'GSE66351':{'col':"Diagnosis",
                                 'ad_val':"diagnosis: AD",
                                 },
                     'GSE134379': {'col':"Diagnosis",
                                 'ad_val':"diagnosis: AD",
                                 },
                     'GSE105109': {'col':"Diagnosis",
                                 'ad_val':"post-mortem diagnosis: Alzheimer's disease",
                                 },
                     'GSE66351_half': {'col': "Diagnosis",
                                  'ad_val': "diagnosis: AD",
                                  },
                     'GSE134379_half': {'col': "Diagnosis",
                                   'ad_val': "diagnosis: AD",
                                   },
                     'GSE105109_half': {'col': "Diagnosis",
                                   'ad_val': "post-mortem diagnosis: Alzheimer's disease",
                                   },
                     }

def read_data(data:Union[str, pd.DataFrame], index_col:Union[int, str, None]=0)->pd.DataFrame:
    if isinstance(data, str):
        if index_col:
            data = pd.read_csv(data, index_col=index_col)
        else:
            data = pd.read_csv(data)
        return data
    elif isinstance(data, pd.DataFrame):
        return data
    else:
        raise TypeError(f'data needs to be a path (str) or a pd.DataFrame {type(data)} was provided')

def createDataSplits(meth_path:Union[str, pd.DataFrame],
                     umeth_path:Union[str, pd.DataFrame],
                     beta_path:Union[str, pd.DataFrame],
                     pheno_path:Union[str, pd.DataFrame],
                     output_path:str,
                     identifier:str,
                     small_design_path:Union[str, pd.DataFrame],
                     diagnosis_col:Union[str, None]=None,
                     ad_diagnosis:Union[str, None]=None,
                     full_design_path:Union[str, pd.DataFrame, None]=None,
                     save_local:bool=True,
                     small_design_local_path:Union[str, pd.DataFrame, None]=None,
                     full_design_local_path:Union[str, pd.DataFrame, None]=None,
                     distortion:str='balanced'):
    # read in the data
    meth = read_data(meth_path, index_col='Unnamed: 0')
    umeth = read_data(umeth_path, index_col='Unnamed: 0')
    beta = read_data(beta_path, index_col='Unnamed: 0')
    small_design = read_data(small_design_path, index_col= "Sample_ID")

    if full_design_local_path:
        full_design = read_data(full_design_path)
    else:
        full_design = None
    pheno = read_data(pheno_path, index_col="Sample_ID")
    if not ad_diagnosis:
        if identifier in label_information:
            ad_diagnosis = label_information[identifier]['ad_val']
        else:
            raise ValueError('No ad-diagnosis was provided')
    if not diagnosis_col and identifier in label_information:
        if identifier in label_information:
            diagnosis_col = label_information[identifier]['col']
        else:
            raise ValueError('No diagnosis was provided')

    if save_local:
        if small_design_local_path:
            small_design_local = read_data(small_design_local_path, index_col="Sample_ID")
        else:
            raise ValueError('Save local was selected and no local small-design matrix was provided')
        if full_design_local_path:
            full_design_local = read_data(full_design_local_path, index_col="Sample_ID")
        else:
            full_design_local = None
    else:
        small_design_local = None
        full_design_local = None


    # create the splits
    splits_pheno = pheno.copy(deep=True)
    N_total = splits_pheno.shape[0] * 1.0 * 1.0
    N_ad = splits_pheno.loc[splits_pheno[diagnosis_col] == ad_diagnosis, :].shape[0] * 1.0 * 1.0
    random_state = 42
    seed(random_state)
    n_splits = 3
    if distortion_ratios[identifier][distortion]['sizes']:
        sizes = distortion_ratios[identifier][distortion]['sizes']
    else:
        raise AttributeError(f'sizes not found in distortion ratios for {identifier} with distortion {distortion}')
    if distortion_ratios[identifier][distortion]['ad_freqs']:
        ad_freqs = distortion_ratios[identifier][distortion]['ad_freqs']
    else:
        raise AttributeError(f'sizes not found in distortion ratios for {identifier} with distortion {distortion}')

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
    ad = set(pheno.loc[pheno[diagnosis_col] == ad_diagnosis, :].index.values)
    print(len(ad))
    other = set(pheno.index.values).difference(ad)  # .difference(fem)
    for i in range(0, n_splits - 1):
        b = set(sample(ad, n_ad[i]))
        ad = ad.difference(b)
        o = set(sample(other, Sizes[i] - n_ad[i]))
        other = other.difference(o)
        sele_samples = b | o
        splits_pheno.loc[sele_samples, "split"] = "Split_" + str(i + 1)
        splits_pheno["Split_" + str(i + 1)] = 0
        splits_pheno.loc[sele_samples, "Split_" + str(i + 1)] = 1

    output_dir = os.path.join(output_path, f"{identifier}_{distortion}_splits")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for i in range(n_splits):
        s = "Split_" + str(i + 1)
        samples = sorted(splits_pheno.loc[splits_pheno["split"] == s, :].index.values)
        print(len(samples))
        splits_pheno.loc[samples, :].to_csv(output_dir + "/" + s + "_pheno.csv")
        meth.loc[:, samples].to_csv(output_dir + "/" + s + "_methylated.csv")
        umeth.loc[:, samples].to_csv(output_dir + "/" + s + "_unmethylated.csv")
        beta.loc[:, samples].to_csv(output_dir + "/" + s + "_betas.csv")
        small_design.loc[samples, :].to_csv(output_dir + "/" + s + "_design.csv")
        if full_design:
            full_design.loc[samples, :].to_csv(output_dir + "/" + s + "_Full_design.csv")
        if save_local:
            small_design_local.loc[samples, :].to_csv(output_dir + "/" + s + "_design_local.csv")
            if full_design_local:
                full_design_local.loc[samples, :].to_csv(output_dir + "/" + s + "_Full_design_local.csv")

    if save_local:
        if small_design_local_path:
            central_design = small_design_local.copy()
            central_design["Cohort_effect"] = splits_pheno["split"]
            central_design["Split1"] = 0
            central_design.loc[central_design["Cohort_effect"] == "Split_1", "Split1"] = 1
            central_design["Split2"] = 0
            central_design.loc[central_design["Cohort_effect"] == "Split_2", "Split2"] = 1
            central_design.drop(columns="Cohort_effect", inplace=True)
            central_design.to_csv(os.path.join(output_dir + "/" + "central_design_matrix.csv"))
        if full_design_local_path:
            central_full_design = full_design_local.copy()
            central_full_design["Cohort_effect"] = splits_pheno["split"]
            central_full_design["Split1"] = 0
            central_full_design.loc[central_full_design["Cohort_effect"] == "Split_1", "Split1"] = 1
            central_full_design["Split2"] = 0
            central_full_design.loc[central_full_design["Cohort_effect"] == "Split_2", "Split2"] = 1
            central_full_design.drop(columns="Cohort_effect", inplace=True)
            central_full_design.to_csv(os.path.join(output_dir + "/" + "full_central_design_matrix.csv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--meth_path",)
    parser.add_argument("-u", "--umeth_path",)
    parser.add_argument("-b", "--beta_path",)
    parser.add_argument("-i", "--identifier",)
    parser.add_argument("-s", "--small_design_path",)
    parser.add_argument("-f", "--full_design_path", default=None)
    parser.add_argument("-d", "--distortion",)
    parser.add_argument("-o", "--output_path",)
    parser.add_argument("-p", "--pheno_path",)
    parser.add_argument('--diagnosis_col')
    parser.add_argument('--ad_diagnosis')
    parser.add_argument("--save_local",action="store_true", default=False)
    parser.add_argument("--small_design_local_path", default=None,)
    parser.add_argument("--full_design_local_path", default=None,)
    args = parser.parse_args()
    if not args.save_local:
        createDataSplits(meth_path=args.meth_path,
                         umeth_path=args.umeth_path,
                         beta_path=args.beta_path,
                         output_path=args.output_path,
                         distortion=args.distortion,
                         pheno_path=args.pheno_path,
                         diagnosis_col=args.diagnosis_col,
                         ad_diagnosis=args.ad_diagnosis,
                         save_local=args.save_local,
                         identifier=args.identifier,
                         small_design_path=args.small_design_path,
                         full_design_path=args.full_design_path)
    elif args.save_local:
        createDataSplits(meth_path=args.meth_path,
                         umeth_path=args.umeth_path,
                         beta_path=args.beta_path,
                         output_path=args.output_path,
                         distortion=args.distortion,
                         pheno_path=args.pheno_path,
                         diagnosis_col=args.diagnosis_col,
                         ad_diagnosis=args.ad_diagnosis,
                         save_local=args.save_local,
                         identifier=args.identifier,
                         small_design_path=args.small_design_path,
                         full_design_path=args.full_design_path,
                         small_design_local_path=args.small_design_local_path,
                         full_design_local_path=args.full_design_local_path)
