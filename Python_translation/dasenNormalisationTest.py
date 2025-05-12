#%%
import pandas as pd
import numpy as np
import sys
try:
    from Federated_Differential_Methylation_Analysis.Python_translation import dasen_normalisation
except ModuleNotFoundError:
    sys.path.append("/cosybio/project/vanElferen/FedEWAS")
    from Federated_Differential_Methylation_Analysis.Python_translation import dasen_normalisation
#%%
np.seterr(under='ignore')
#%%
# read in the data
unmeth = pd.read_csv('/cosybio/project/vanElferen/FedEWAS/Data/QC_GSE66351/Filtered_Unmethylated.csv',
                     index_col=0)
meth = pd.read_csv('/cosybio/project/vanElferen/FedEWAS/Data/QC_GSE66351/Filtered_Methylated.csv',
                   index_col=0)
annotation = pd.read_csv('/cosybio/project/vanElferen/FedEWAS/Data/GSE134379_Raw/GPL13534_HumanMethylation450_15017482_v.1.1.csv',
                         skiprows=7, index_col=0, low_memory=False)
annotation_data = annotation.loc[list(set(meth.index.values).intersection(set(annotation.index.values))), :]
#%%
# run the complete normalisation thing
normalised_betas = dasen_normalisation.dasen_normalisation(unmeth, meth, annotation_data.loc[:,"Infinium_Design_Type"])
#%%
# save the normalised beta for comparison again R dasen implementation results
normalised_betas.to_csv('/cosybio/project/vanElferen/FedEWAS/Data/QC_Python/new_dfsfit_implementation_normalised_betas.csv')
#%%