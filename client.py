from sys import stderr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.distributions.mixture_rvs import mixture_rvs
import statsmodels.api as sm
import scipy
import scipy.stats
from scipy.stats import rankdata
from scipy import interpolate
import re
import sys

# dasen normalisation local substeps
def dfs2_python(x, probe_type):
    import statsmodels.api as sm
    from statsmodels.distributions.mixture_rvs import mixture_rvs

    # new code version that should work on one column at a time
    x_copy = x.copy()
    KD_one = sm.nonparametric.KDEUnivariate(x_copy[probe_type.squeeze() == "I"])
    KD_one.fit(gridsize=5000)
    one = int(KD_one.support[np.where(np.max(KD_one.density))])
    KD_two = sm.nonparametric.KDEUnivariate(x_copy[probe_type.squeeze() == "II"])
    KD_two.fit(gridsize=5000)
    two = int(KD_two.support[np.where(np.max(KD_two.density))])
    out = np.max(one) - np.max(two) #not quite sure if any of this is correct
    return out

def dfsfit_python(x, probe_type):
    import statsmodels.api as sm
    import re
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
    
    fit_dist = sm.OLS.from_formula("dis_diff ~ scol + srow", dis_diff).fit()
    dis_diff = [fit_dist.fittedvalues]
    n = probe_type.squeeze() == "I"
    tI_correction = np.tile(np.array(dis_diff), (sum(n),1))
    x[probe_type.squeeze() == "I"] = x[probe_type.squeeze() == "I"] - tI_correction
    return x

class Client:
    def __init__(self, cohort_name, design_matrix_filepath, methylated_filepath, unmethylated_filepath, probe_annotation_path):
        self.cohort_name = cohort_name
        self.probes = None
        self.sample_names = None
        self.read_data(design_matrix_filepath, methylated_filepath, unmethylated_filepath, probe_annotation_path)
        
        #data slots
        self.raw_methylated = None
        self.raw_unmethylated = None
        self.probe_annotation = None
        self.designmatrix = None
        self.designcolumns = None

        # (dasen) normalisation
        self.methylated_dist = None
        self.unmethylated_dist = None
        self.methI_local_norm_par = None
        self.methII_local_norm_par = None
        self.unmethI_local_norm_par = None
        self.unmethII_local_norm_par = None
        self.methnorm = None
        self.unmethnorm = None
        self.betas = None

        # EWAS
        self.xtx = None
        self.xty = None
        self.coef = None
        self.stnd_err = None
        self.p_value = None
        self.EWAS_results = None

    def read_data(self, design_matrix_filepath, methylated_filepath, unmethylated_filepath, probe_annotation_path):
        '''
        Read in the data provided to the client which will contain:
            The design matrix: numerical variable, binary representation of two-level categorical variable
            or binary dummy representation of multi-level categorical variable
        Checks that the samples included in the design matrix are also present in the data and vice-versa and creates 
        a design matrix containg only samples that are present in the data and removes any samples that are not present
        in the design matrix from the data
        '''
        # data
        self.raw_methylated = pd.read_csv(methylated_filepath, index_col=0)
        self.raw_unmethylated = pd.read_csv(unmethylated_filepath, index_col=0)
        self.probe_annotation = pd.read_csv(probe_annotation_path, index_col=0)

        # design
        self.designmatrix = pd.read_csv(design_matrix_filepath, index_col=0)
        self.designcolumns = list(self.designmatrix.columns.values)

        # check that the indexes of the methylated and unmethylated dataframes are the same
        # and use one of them to set the probes attribute
        if np.all(self.raw_methylated.index == self.raw_unmethylated.index):
            self.probes = list(self.raw_methylated.index.values)
            
        else:
            print("Probe names in the methylated and unmethylated files don't match", file=sys.stderr)
            exit(1)

        # set the sample names - sanity check that the same samples are in the methylated and unmethylated files
        if np.all(self.raw_methylated.columns.values == self.raw_unmethylated.columns.values):
            self.sample_names = list(self.raw_methylated.columns.values)
        
        # make sure that the samples appear in the same order in the data matrix as in the design matrix
        design_rows = set(self.designmatrix.index.values)
        data_columns = set(self.raw_methylated.columns.values)
        if not np.all(design_rows == data_columns):
            print("The rows of the design matrix and the columns of the data files don't match", file=stderr)
            # keep the samples that are present in both the data and the design matrix
            samples_keep = sorted(design_rows.intersection(data_columns))
            samples_not_in_design = design_rows.difference(samples_keep)
            print("These samples are not present in the design matrix but are present in the data: %s"%(samples_not_in_design), file=stderr)
            samples_not_in_data = data_columns.difference(samples_keep)
            print("These samples are not present in the data but are present in the design matrix: %s"%(samples_not_in_data), file=stderr)
            
            #select the overlapping samples from the data and the design matrix and save them in the correct attribute slots
            self.designmatrix = self.designmatrix.loc[samples_keep,:]
            self.raw_methylated = self.raw_methylated.loc[:,samples_keep]
            self.raw_unmethylated = self.raw_unmethylated.loc[:,samples_keep]
    
    def input_validation(self, global_conditions):
        ''' Checks the input for the following:
                -if the same probes are included in both the methylated and unmethylated dataframe
                -if the probe type annotation data of all probes in the data is present and remove the annotation
                for any probes not present in the data
                -if all the global conditions (specified at server level) are present in the local design matrix
                -if there are any additional conditions to the globally specified conditions in the local design
                matrix and removes these
        '''      
        # check probe type annotation data
        if not np.all(self.probes == self.probe_annotation.index.values):
            data_probes = set(self.probes)
            annotation_probes = set(self.probe_annotation.index.values)
            no_annotation = data_probes.difference(annotation_probes)
            if no_annotation > 0:
                print("Probe type annotation for %s probes in data is missing"%(len(no_annotation)), file = stderr)
                exit(1)
            no_data = annotation_probes.difference(data_probes)
            if no_data > 0:
                keep_annotation = annotation_probes.intersection(data_probes)
                self.probe_annotation = self.probe_annotation.loc[keep_annotation,:]

        #check the design matrix columns against the globally specified conditions
        local_conditions = set(self.designcolumns)
        global_conditions = set(global_conditions)
        if not np.all(self.designcolumns == global_conditions):
            no_local = global_conditions.difference(local_conditions)
            no_global = local_conditions.difference(global_conditions)
            if len(no_local) > 0:
                print("%s global conditions are not present in the local design matrix"%(len(no_local)), file=stderr)
                exit(1)
            if len(no_global) > 0:
                print("%s conditions in the design matrix are not specified in the global conditions and will be removed"%(len(no_global)), file=stderr)
                conditions_to_keep = list(local_conditions.intersection(global_conditions))
                self.designmatrix = self.designmatrix.loc[:,conditions_to_keep]
                self.designcolumns = list(self.designmatrix.columns.values)
    
    def cohort_effects(self, cohort_names):
        for cohort in cohort_names:
            if self.cohort_name == cohort:
                self.designmatrix[cohort] = 1
            else:
                self.designmatrix[cohort] = 0
        self.designcolumns = list(self.designmatrix.columns.values)
        
    # client level computation for (dasen)normalisation
    def intensity_distributions(self):
        self.methylated_dist = dfsfit_python(self.raw_methylated, self.probe_annotation)
        self.unmethylated_dist = dfsfit_python(self.raw_unmethylated, self.probe_annotation)
        return self.methylated_dist, self.unmethylated_dist
    
    def local_normalisation_parameters(self):
        # for the methylated type I probes, methylated type II, unmethylated type I and unmethylated type II        
        methI_data = self.methylated_dist[self.probe_annotation.values == "I"].to_numpy() # not working as expected
        methII_data = self.methylated_dist[self.probe_annotation.values == "II"].to_numpy()
        unmethI_data = self.unmethylated_dist[self.probe_annotation.values == "I"].to_numpy()
        unmethII_data = self.unmethylated_dist[self.probe_annotation.values == "II"].to_numpy()
        data_list = [methI_data, methII_data, unmethI_data, unmethII_data]
        observations = []
        row_sums = []
        n_columns = []
        for data in data_list:
            n,m = data.shape
            n_observations = np.array(m * [n]) # create an array (1,m) that has the number of observations (n) for each column 

            idx_array = np.empty((n,m))
            idx_array[:] = np.nan
            sorted = idx_array.copy()
            for col in range(0,m):
                col_s = np.sort(data[:,col]) # sort the column data
                #create a dataframe with the indices
                col_id = data[:,col]
                bla = np.argsort(col_id)
                col_idx = np.empty_like(col_id)
                col_idx[bla] = np.arange(len(col_id))

                obs_col = len(col_s) - np.count_nonzero(np.isnan(data))
                if obs_col < n:
                    None
                else:
                    sorted[:,col] = col_s
                    idx_array[:,col] = col_idx
            row_sum = np.sum(sorted, axis = 1)
            observations.append(n_observations)
            row_sums.append(row_sum)
            n_columns.append(m)
        self.methI_local_norm_par = [observations[0], row_sums[0], n_columns[0]]
        self.methII_local_norm_par = [observations[1], row_sums[1], n_columns[1]]
        self.unmethI_local_norm_par = [observations[2], row_sums[2], n_columns[2]]
        self.unmethII_local_norm_par = [observations[3], row_sums[3], n_columns[3]]

        return self.methI_local_norm_par, self.methII_local_norm_par, self.unmethI_local_norm_par, self.unmethII_local_norm_par
        #return observations, row_sums, n_columns
        

    def final_normalisation(self, probe_type_means):
        methI_data = self.methylated_dist[self.probe_annotation.values == "I"].to_numpy()
        methII_data = self.methylated_dist[self.probe_annotation.values == "II"].to_numpy()
        unmethI_data = self.unmethylated_dist[self.probe_annotation.values == "I"].to_numpy()
        unmethII_data = self.unmethylated_dist[self.probe_annotation.values == "II"].to_numpy()
        data_list = [methI_data, methII_data, unmethI_data, unmethII_data]

        data_out = []
        for i in range(0,len(data_list)):
            n,m = data_list[i].shape
            k = np.arange(n)/(n-1)
            results = np.empty((n,m))
            results[:] = np.nan
            for col in range(0,m):
                rank = rankdata(data_list[i][:,col], method="average")
                f = scipy.interpolate.interp1d(k, probe_type_means[i])
                results[:,col] = f((rank - 1)/(n-1))
            data_out.append(pd.DataFrame(results))
        #save the normalised methylated intensities
        data_out[0].set_index(self.methylated_dist[self.probe_annotation.squeeze() =="I"].index, inplace = True)
        data_out[0].columns = self.sample_names

        data_out[1].set_index(self.methylated_dist[self.probe_annotation.squeeze() =="II"].index, inplace = True)
        data_out[1].columns = self.sample_names
        self.methnorm = pd.concat([data_out[0], data_out[1]])
        self.methnorm.sort_index(inplace=True)

        #save the normalised unmethylated intensities
        data_out[2].set_index(self.methylated_dist[self.probe_annotation.squeeze() =="I"].index, inplace = True)
        data_out[2].columns = self.sample_names

        data_out[3].set_index(self.methylated_dist[self.probe_annotation.squeeze() =="II"].index, inplace = True)
        data_out[3].columns = self.sample_names
        self.unmethnorm = pd.concat([data_out[2], data_out[3]])
        self.unmethnorm.sort_index(inplace=True)
        self.betas = self.methnorm/(self.methnorm + self.unmethnorm + 100)
        return self.betas

    # client level computations for EWAS based on simple linear model
    def local_xtx_xty(self):
        '''
        Calculate the local intermediate matrices that are sent to the server
        where they are used to calculate the global intermediates
        '''
        x_matrix = self.designmatrix.values
        y_matrix = self.betas.values
        n = y_matrix.shape[0] # number of genes
        m = x_matrix.shape[1] # number of conditions
        self.xtx = np.zeros((n,m,m)) # create zeroes array with space to hold the xt_x matrix (m,m) for each probe (n)
        self.xty = np.zeros((n,m)) # create zeroes array with space to hold the xt_y matrix (1,m) for each probe (n)
        for i in range(0,n):
            y = y_matrix[i,:]
            self.xtx[i,:,:] = x_matrix.T @ x_matrix
            self.xty[i,:] = x_matrix.T @ y
        return self.xtx, self.xty
    
    def calculate_EWAS_results(self, global_xtx, global_xty):
        n,m = global_xty.shape
        self.coef = np.zeros((n,m))
        self.stnd_err = np.zeros((n,m))
        self.p_value = np.zeros((n,m))

        for i in range(0,n):
            xtx_inv = np.linalg.inv(global_xtx[i,:,:])
            self.coef[i,:] = xtx_inv @ global_xty[i,:]
            self.stnd_err[i,:] = np.diag(xtx_inv)
            t = self.coef[i,:]/self.stnd_err[i,:]
            df = m-2
            self.p_value[i,:] = scipy.stats.t.sf(t, df)
        # create EWAS results dataframe with all information grouped by variable/confounder
        coef = pd.DataFrame(self.coef,index=self.probes, columns= self.designcolumns)
        stnErr = pd.DataFrame(self.stnd_err, index=self.probes, columns= self.designcolumns)
        p_val = pd.DataFrame(self.p_value, index=self.probes, columns= self.designcolumns)
        self.EWAS_results = pd.concat([coef, stnErr, p_val], axis = 1, keys = ["Coefficient", "StandardError", "P-value"])
        
        return self.EWAS_results


        

