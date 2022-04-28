from threading import local
import pandas as pd
import numpy as np
#import scipy
#from scipy.stats import rankdata
#from scipy import interpolate
import re
import sys

class Server:
    def __init__(self, variables, confounders):
        self.variables = variables
        self.confounders = confounders
        self.client_names = []
        self.samples_per_client = []
        self.global_probes = []

        self.global_xtx = None
        
    # get client information
    def get_clients(self, cohort_name, client_probes, client_n_samples, client_design_columns):
        
        None

    # aggragate local means for quantile normalisation
    def aggregate_QN_means(*local_sums):
        probe_type_means = []
        for i in range(0,len(local_sums)):
            n = len(local_sums[i])
            if n == 0:
                print("Requires at least 1 local sum list to function", file=sys.stderr)
            elif n == 1:
                global_means = input[0][1] # return row means back to the client
            else:
                s = local_sums[0][1]
                obs = local_sums[0][2]
                for j in range(1, len(local_sums)):
                    s = s + local_sums[i][1]
                    obs = obs + local_sums[i][2]
                    global_means = obs/s
            probe_type_means.append(global_means)
        return probe_type_means
    
    # aggregate local xt_x to global xt_x matrix for linear regression
    def global_regression_parameter(self, *local_xtx):
        n = self.global_probes
        m = (self.variables + self.confounders)
        self.global_xtx = np.zeros((n,m,m))
        for matrix in local_xtx:
            self.global_xtx += local_xtx[matrix]
        return self.global_xtx

