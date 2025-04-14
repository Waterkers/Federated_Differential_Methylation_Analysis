import pandas as pd
from random import seed, sample
import os
import argparse

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
                     }

def createDataSplits(meth_path:str,
                     umeth_path:str,
                     beta_path:str,
                     output_path:str,
                     identifier:str,
                     small_design_path:str,
                     full_design_path:str,
                     pheno_path:str,
                     diagnosis_col:str,
                     ad_diagnosis:str,
                     save_local:bool=True,
                     small_design_local_path:str=None,
                     full_design_local_path:str=None,
                     distortion:str='balanced'):
    meth = pd.read_csv(meth_path)
    umeth = pd.read_csv(umeth_path)
    beta = pd.read_csv(beta_path)
    small_design = pd.read_csv(small_design_path)
    full_design = pd.read_csv(full_design_path)
    pheno = pd.read_csv(pheno_path)
    splits_pheno = pheno.copy(deep=True)
    if save_local:
        if small_design_local_path:
            small_design_local = pd.read_csv(small_design_local_path)
        else:
            small_design_local = None
        if full_design_local_path:
            full_design_local = pd.read_csv(full_design_local_path)
        else:
            full_design_local = None

    N_total = splits_pheno.shape[0] * 1.0 * 1.0
    N_ad = splits_pheno.loc[splits_pheno[diagnosis_col] == ad_diagnosis, :].shape[0] * 1.0 * 1.0
    random_state = 42
    seed(random_state)
    n_splits = 3
    if distortion_ratios[identifier][distortion]['sizes'] & distortion_ratios[identifier][distortion]['ad_freqs']:
        sizes = distortion_ratios[identifier][distortion]['sizes']
        ad_freqs = distortion_ratios[identifier][distortion]['ad_freqs']
    else:
        raise AttributeError(f'sizes or ad_freqs not found in distortion ratios for {identifier} with distortion {distortion}')

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

    output_dir = os.path.join(output_path, f"{identifier}_strong_splits")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for i in range(n_splits):
        s = "Split_" + str(i + 1)
        samples = sorted(splits_pheno.loc[splits_pheno["split"] == s, :].index.values)
        splits_pheno.loc[samples, :].to_csv(output_dir + "/" + s + "_pheno.csv")
        meth.loc[:, samples].to_csv(output_dir + "/" + s + "_methylated.csv")
        umeth.loc[:, samples].to_csv(output_dir + "/" + s + "_unmethylated.csv")
        beta.loc[:, samples].to_csv(output_dir + "/" + s + "_betas.csv")
        small_design.loc[samples, :].to_csv(output_dir + "/" + s + "_design.csv")
        full_design.loc[samples, :].to_csv(output_dir + "/" + s + "_Full_design.csv")
        if save_local:
            small_design_local.loc[samples, :].to_csv(output_dir + "/" + s + "_design_local.csv")
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
    parser.add_argument("-f", "--full_design_path",)
    parser.add_argument("-d", "--distortion",)
    parser.add_argument("-o", "--output_path",)
    parser.add_argument("-p", "--pheno_path",)
    parser.add_argument('--diagnosis_col', default='Diagnosis')
    parser.add_argument('--ad_diagnosis', default=' AD')
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
        if args.small_design_local_path & args.full_design_local_path:
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
        else:
            raise ValueError('Selected save local but not path to local design matrices provided')