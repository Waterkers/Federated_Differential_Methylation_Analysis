from sys import stderr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.distributions.mixture_rvs import mixture_rvs
from statsmodels.stats.multitest import multipletests
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
    from statsmodels.nonparametric.kde import KDEUnivariate
    import DensityR
    from scipy.stats import scoreatpercentile


    # new code version that should work on one column at a time
    x_copy = x.copy()
    # KDEUnivariate attempt to match r
    KD_one = DensityR.KDEUnivariate_rDensity(x_copy[probe_type.squeeze() == "I"])
    KD_one.fit(gridsize=2**15, low=0, high=5000)
    one = KD_one.support[np.argmax(KD_one.density)]

    # KDEUnivariate attempt to match r
    KD_two = DensityR.KDEUnivariate_rDensity(x_copy[probe_type.squeeze() == "II"])
    KD_two.fit(gridsize=2**15, low=0, high=5000)
    two = KD_two.support[np.argmax(KD_two.density)]
    
    out = np.max(one) - np.max(two) #not quite sure if any of this is correct
    return out

def dfsfit_python(x, probe_type):
    import statsmodels.api as sm
    import re
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
    x[probe_type.squeeze() == "I"] = x[probe_type.squeeze() == "I"] - tI_correction
    return x

class Client:
    def __init__(self, cohort_name, design_matrix_filepath, methylated_filepath, unmethylated_filepath, probe_annotation_path):
        self.cohort_name = cohort_name
        self.probes = None
        self.sample_names = None
                
        #data slots
        self.raw_methylated = None
        self.raw_unmethylated = None
        self.probe_annotation = None
        self.designmatrix = None
        self.designcolumns = None
        self.read_data(design_matrix_filepath, methylated_filepath, unmethylated_filepath, probe_annotation_path)

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
        self.mu = None
        self.xtx = None
        self.xty = None
        self.SSE = None
        self.cov_coef = None
        """ self.p_value = None
        self.EWAS_results = None
        self.coef = None
        self.stnd_err = None """

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
            sys.exit(1)

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
    
    def input_validation(self, global_conditions, global_probes):
        ''' Checks the input for the following:
                -if the probe type annotation data of all probes in the data is present and remove the annotation
                for any probes not present in the data
                - retain only the probes that are present in all datasets analysed in the project
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
                sys.exit(1)
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
                sys.exit(1)
            if len(no_global) > 0:
                print("%s conditions in the design matrix are not specified in the global conditions and will be removed"%(len(no_global)), file=stderr)
                conditions_to_keep = list(local_conditions.intersection(global_conditions))
                self.designmatrix = self.designmatrix.loc[:,conditions_to_keep]
                self.designcolumns = list(self.designmatrix.columns.values)
        
        #check the local probes against the global probes and remove any local probes that are not present in the global dataset of the project
            # could later be changed to impute the missing probes in each dataset based on the values of these in the datasets that do have them?
        local_probes = set(self.probes)
        global_probes = set(global_probes)
        if not np.all(local_probes == global_probes):
            not_local = global_probes.difference(local_probes)
            not_global = local_probes.difference(global_probes)
            if len(not_local) > 0:
                print("%s of global probes not present in local dataset - check global probe definition"%(len(not_local)), file=stderr)
            if len(not_global) > 0:
                print("%s probes are not present in the global dataset and will be removed before analysis"%(len(not_global)), file=stderr)
                probes_to_keep = list(local_probes.intersection(global_probes))
                self.raw_methylated = self.raw_methylated.loc[probes_to_keep, :]
                self.raw_unmethylated = self.raw_unmethylated.loc[probes_to_keep, :]
                self.probes = probes_to_keep

    def find_unique_SentrixIDS(self):
        sentrix_ids = self.designmatrix.loc[:, "Sentrix_ID"]
        self.unique_SentrixIDS = list(set(sentrix_ids))
    
    def find_unique_PlateIDS(self):
        plate_ids = self.designmatrix.loc[:, "plate_id"]
        self.unique_PlateIDS = list(set(plate_ids))

    def cohort_effects(self, cohort_effects):
        for cohort in cohort_effects:
            if self.cohort_name == cohort:
                self.designmatrix[cohort] = 1
            else:
                self.designmatrix[cohort] = 0
        self.designcolumns = list(self.designmatrix.columns.values)
        
    def SentrixID_effects(self, global_SentrixID):
        for ID in global_SentrixID:
            #if ID in self.unique_SentrixIDS:
            self.designmatrix[ID] = 0
            self.designmatrix[ID].loc[self.designmatrix["Sentrix_ID"] == ID] = 1
            #self.designmatrix[ID].loc[self.designmatrix["Sentrix_ID"] != ID] = 0
        self.designmatrix.drop(columns="Sentrix_ID", inplace=True)
        self.designcolumns = list(self.designmatrix.columns.values)

    def PlateID_effects(self, global_PlateID):
        for ID in global_PlateID:
            #if ID in self.unique_PlateIDS:
            self.designmatrix[ID] = 0
            self.designmatrix[ID].loc[self.designmatrix["Plate_ID"] == ID] = 1
        self.designmatrix.drop(columns="Plate_ID", inplace=True)
        self.designcolumns = list(self.designmatrix.columns.values)

    # client level computation for (dasen)normalisation
    def intensity_distributions(self):
        self.methylated_dist = dfsfit_python(self.raw_methylated, self.probe_annotation)
        self.unmethylated_dist = dfsfit_python(self.raw_unmethylated, self.probe_annotation)
        return self.methylated_dist, self.unmethylated_dist
    
    def local_normalisation_parameters(self):
        # for the methylated type I probes, methylated type II, unmethylated type I and unmethylated type II        
        methI_data = self.methylated_dist[self.probe_annotation.values == "I"].copy().to_numpy() # not working as expected
        methII_data = self.methylated_dist[self.probe_annotation.values == "II"].copy().to_numpy()
        unmethI_data = self.unmethylated_dist[self.probe_annotation.values == "I"].copy().to_numpy()
        unmethII_data = self.unmethylated_dist[self.probe_annotation.values == "II"].copy().to_numpy()
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

        #return self.methI_local_norm_par, self.methII_local_norm_par, self.unmethI_local_norm_par, self.unmethII_local_norm_par
        #return observations, row_sums, n_columns
        

    def final_normalisation(self, probe_type_means):
        methI_data = self.methylated_dist[self.probe_annotation.values == "I"].copy().to_numpy()
        methII_data = self.methylated_dist[self.probe_annotation.values == "II"].copy().to_numpy()
        unmethI_data = self.unmethylated_dist[self.probe_annotation.values == "I"].copy().to_numpy()
        unmethII_data = self.unmethylated_dist[self.probe_annotation.values == "II"].copy().to_numpy()
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
    def local_xtx_xty(self, weighted = False):
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
        
        if weighted:
            weighted = True
            W = np.sqrt(self.weights)
            y_matrix = np.multiply(y_matrix,W) 

        for i in range(0,n):
            y = y_matrix[i,:]
            if weighted:
                Xw = np.multiply(x_matrix,W[i,:].reshape(-1, 1)) # algebraic multiplications by W
                self.XtX[i,:,:] = Xw.T @ Xw 
                self.XtY[i,:] = Xw.T @ y  
            else: 
                self.xtx[i,:,:] = x_matrix.T @ x_matrix
                self.xty[i,:] = x_matrix.T @ y
        return self.xtx, self.xty
    
    def compute_SSE_and_cov_coef(self,beta, weighted = False):
        x_matrix = self.designmatrix.values
        y_matrix = self.betas.values 
        n = y_matrix.shape[0]
        m = x_matrix.shape[1]
        self.SSE = np.zeros(n)
        self.mu = np.zeros(m)
        if weighted:
            W = np.sqrt(self.weights)
            y_matrix = np.multiply(y_matrix,W)
        for i in range(0,n): 
            y = y_matrix[i,:]
            if weighted:
                Xw = np.multiply(x_matrix,W[i,:].reshape(-1, 1))
                self.mu[i,] =  Xw @ beta[i,:] # fitted logCPM 
            else:
                self.mu[i,] =  x_matrix @ beta[i,:] # fitted logCPM 

            self.SSE[i] = np.sum((y - self.mu[i,])**2) # local SSE
        #print("mu:",self.mu.shape)
        Q,R = np.linalg.qr(x_matrix)
        self.cov_coef = R.T @ R
        self.cov_coef = x_matrix.T @ x_matrix
        return self.SSE, self.cov_coef

    """ def calculate_EWAS_results(self, global_xtx, global_xty):
        n,m = global_xty.shape
        self.coef = np.zeros((n,m))
        self.stnd_err = np.zeros((n,m))
        self.p_value = np.zeros((n,m))
        self.p_value_cor = np.zeros((n,m))
        self.mchange = np.zeros((n,m))
        
        for i in range(0,n):
            xtx_inv = np.linalg.inv(global_xtx[i,:,:])
            self.coef[i,:] = xtx_inv @ global_xty[i,:]
            self.stnd_err[i,:] = np.diag(xtx_inv)
            t = self.coef[i,:]/self.stnd_err[i,:]
            df = m-2
            self.p_value[i,:] = scipy.stats.t.sf(t, df)
            #self.p_value_cor[i,:] = multipletests(self.p_value[i,:], method="fdr_bh")[1]
        # create EWAS results dataframe with all information grouped by variable/confounder
        for i in range(0,m):
            self.p_value_cor[:,i] = multipletests(self.p_value[:,i], method="fdr_bh")[1]
        
        for i in range(0,m):
            self.mchange[:,i] = self.coef[:,i]
        coef = pd.DataFrame(self.coef,index=self.probes, columns= self.designcolumns)
        stnErr = pd.DataFrame(self.stnd_err, index=self.probes, columns= self.designcolumns)
        p_val = pd.DataFrame(self.p_value, index=self.probes, columns= self.designcolumns)
        p_val_corrected = pd.DataFrame(self.p_value_cor, index=self.probes, columns= self.designcolumns)
        mchange = pd.DataFrame(self.mchange, index = self.probes, columns=self.designcolumns)
        # create a dataframe with the corrected p-values
        
        self.EWAS_results = pd.DataFrame([coef["Diagnosis"], stnErr["Diagnosis"], p_val["Diagnosis"], p_val_corrected["Diagnosis"], mchange["Diagnosis"]],
            columns=["Coefficient", "Standar Error", "P-value", "Corrected P-value", "Methylation Change"])
        
        return self.EWAS_results """


        


