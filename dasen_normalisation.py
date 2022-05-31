import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.distributions.mixture_rvs import mixture_rvs
import DensityR
import statsmodels.api as sm
import re

def dasen_normalisation(unmethylated, methylated, probe_type, base = 100):
    """
    computes the dasen normalised beta values: quantile normalises the unmethylated and methylated intensities, per probe type,
    and uses these normalised intensities to calculate the beta values

    Input arguments:
    unmethylated = dataframe of unmethylated intensities
    methylated = dataframe of methylated intensities
    probe_type = series indicating the type of each probe (Type I or Type II)

    Returns: a dataframe of normalised beta values
    """
    # fit the probability ditribution to the methylated and unmethylated probe intensities based on their probe type
    unmethylated_fit = dfsfit_python(unmethylated, probe_type)
    methylated_fit = dfsfit_python(methylated, probe_type)

    # calculate the quantile normalised values for the methylated and unmethylated probe intensities based on the estimated distribution values
    unmethylated.loc[probe_type == "I"] = quantile_normalise(unmethylated_fit.loc[probe_type == "I"])
    unmethylated.loc[probe_type == "II"] = quantile_normalise(unmethylated_fit.loc[probe_type == "II"])

    methylated.loc[probe_type == "I"] = quantile_normalise(methylated_fit.loc[probe_type == "I"])
    methylated.loc[probe_type == "II"] = quantile_normalise(methylated_fit.loc[probe_type == "II"])

    # calculate the new beta values based on the per probe normalised methylated and unmethylated probe intentisity values
    betas = methylated/(methylated + unmethylated + base) 
    return betas


def dfs2_python(x, probe_type):
    
    # new code version that should work on one column at a time
    x_copy = x.copy()
    KD_one = DensityR.KDEUnivariate_rDensity(x_copy[probe_type == "I"])
    KD_one.fit(gridsize=2**15, low=0, high=5000)
    one = int(KD_one.support[np.argmax(KD_one.density)])

    KD_two = DensityR.KDEUnivariate_rDensity(x_copy[probe_type == "II"])
    KD_two.fit(gridsize=2**15, low=0, high=5000)
    two = int(KD_two.support[np.argmax(KD_two.density)])

    out = np.max(one) - np.max(two) 
    return out

def dfsfit_python(x, probe_type):
    
    x = x.copy()
    dis_diff = x.apply(dfs2_python, args = (probe_type,), axis=0) #create a dataframe/array of the values when dfs2 is applied to each column
    
    roco = []
    for col_name in x.columns.values.tolist() :
        found = re.search("(R0[1-9]C0[1-9])", col_name).group(1)
        roco.append(found) 
    
    srow = []
    scol = []
    for ro in roco:
        row = int(ro[2])
        srow.append(row)
        col = int(ro[5])
        scol.append(col)
    
    roco_zip = list(zip(srow, scol))
    data = pd.DataFrame(roco_zip, index = x.columns.values, columns = ["srow", "scol"])
    data.insert(loc = 0, column="dis_diff", value=dis_diff)

    fit_dist = sm.OLS.from_formula("dis_diff ~ scol + srow", data).fit()
    dis_diff = [fit_dist.fittedvalues]
    n = probe_type.squeeze() == "I"
    tI_correction = np.tile(np.array(dis_diff), (sum(n),1))
    x.loc[probe_type == "I"] = x.loc[probe_type == "I"] - tI_correction
    return x

def quantile_normalise(input_data):
    """
    input_data = a dataframe that needs to be quantile normalised
    returns a quantile normalised version of the input_data
    """
    from scipy import interpolate
    from scipy.stats import rankdata
    data = input_data.copy().to_numpy()
    n,m = data.shape
    sorted = np.empty((n,m))
    sorted[:] = np.nan
    for col in range(0,m):
        col_s = np.sort(data[:,col])
        sorted[:,col] = col_s
    row_means = np.mean(sorted, axis = 1)
    results = np.empty((n,m))
    results[:] = np.nan
    k = np.arange(n)/(n-1)
    for col in range(0,m):
        rank = rankdata(data[:,col], method="average")
        f = interpolate.interp1d(k, row_means)
        results[:,col] = f((rank - 1)/(n-1))
    QN_out = pd.DataFrame(results, index=input_data.index.values, columns=input_data.columns.values)
    """ data_sorted = pd.DataFrame(np.sort(input_data.values, axis = 0), index = input_data.index, columns = input_data.columns) #sort the values of each column (sample) and keep the original row 
    # and column names
    data_sorted_means = data_sorted.mean(axis = 1) # calulate the row means of the sorted data -> these means will be used to replace the raw values in the data
    n, m = input_data.shape
    k = np.arange(n)/(n-1)
    data_sorted_means_int = data_sorted_means.apply(interpolate.interp1d, axis = 0, x=k)
    data_sorted_means.index = np.arange(1, len(data_sorted_means)+1) # this sets the index so it will correspond to the descending ranks that will be assigned to the original 
    # data in the dataframe. This way the row means, which are sorted loweste to highest, can be used to replace the raw data in the correct order
    data_rank = input_data.rank(method = "min").stack().astype(int) # get the rank of the values for each sample in the raw dataset in integer format and change the dataframe so that
    # the columns become the rows, with a multi-index indicating probe as the highest level and the samples for that probe as the second level
    QN_data = data_rank.map(data_sorted_means_int).unstack() # map the row mean values onto the matching ranks obtained from the original dataframe and bring it back to a row = probe
    # and column = sample format """
    return QN_out
