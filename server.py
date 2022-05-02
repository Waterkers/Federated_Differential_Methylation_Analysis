from threading import local
import pandas as pd
import numpy as np
#import scipy
#from scipy.stats import rankdata
#from scipy import interpolate
import re
import sys

class Server:
    def __init__(self, variables:list, confounders:list):
        self.variables = variables
        self.confounders = confounders
        self.client_names = []
        self.samples_per_client = []
        self.client_total_probes = []
        self.global_probes = []

        self.global_xtx = None
        self.global_xty = None
        
    
    def get_clients(self, cohort_name, client_probes, client_n_samples):
        '''
        get client information
        '''
        self.client_names.append(cohort_name)
        self.samples_per_client.append(len(client_n_samples))
        self.client_total_probes.extend(client_probes)
    
           
    def find_global_probes(self):
        '''
        find the probes that are present in the data held by the clients 
        '''
        n_clients = len(self.client_names)
        #common_index = []
        for probe in self.client_total_probes:
            if (self.client_total_probes.count(probe) == n_clients):
                self.global_probes.append(probe)
            else: 
                pass  
        return self.global_probes

    def return_global_conditions(self):
        self.confounders.extend(self.variables)
        global_conditions = self.confounders
        return global_conditions
    
    def aggregate_QN_means(self, *local_sums):
        '''
        aggragate local means for quantile normalisation
        '''
        n = len(local_sums)
        if n == 0:
            print("Requires at least 1 local sum list to function", file=sys.stderr)
        elif n == 1:
            global_means = local_sums[1] # return row means back to the client
        else:
            methI = local_sums[0][0][1]
            methII = local_sums[0][1][1]
            unmethI =  local_sums[0][2][1]
            unmethII = local_sums[0][3][1]

            obs_mI = local_sums[0][0][2]
            obs_mII = local_sums[0][1][2]
            obs_uI = local_sums[0][2][2]
            obs_uII = local_sums[0][3][2]

            #probe_type_means = []
            for i in range(1, len(local_sums)):
                methI = methI + local_sums[i][0][1]
                methII = methII + local_sums[i][1][1]
                unmethI = unmethI + local_sums[i][2][1]
                unmethII = unmethII + local_sums[i][3][1]
        
        
                obs_mI = obs_mI + local_sums[i][0][2]
                obs_mII = obs_mII + local_sums[i][1][2]
                obs_uI = obs_uI + local_sums[i][2][2]
                obs_uII = obs_uII + local_sums[i][3][2]

                gmethI = methI/obs_mI
                gmethII = methII/obs_mII
                gunmethI = unmethI/obs_uI
                gunmethII = unmethII/obs_uII
        
            #probe_type_means.append(global_means)
        return gmethI, gmethII, gunmethI, gunmethII

        """ n = len(local_sums)
        if n == 0:
            print("Requires at least 1 local sum list to function", file=sys.stderr)
        elif n == 1:
            global_means = local_sums[1] # return row means back to the client
        else:
            probe_type_means = []
            for i in range(len(local_sums)):
                methI, methII, unmethI, unmethII = local_sums[i]
                data = [methI, methII, unmethI, unmethII]
        
                s = data[0][1]
                obs = data[0][2]
                for list in range(1,len(data)):
                    s = s + list[1]
                    obs = obs + list[2]
                    global_means = obs/s
                probe_type_means.append(global_means)
        return probe_type_means """
    
    
    def global_regression_parameter(self, *local_xt_matrices):
        '''
        aggregate local xt_x to global xt_x matrix for linear regression
        '''
        n = self.global_probes
        m = (self.variables + self.confounders)
        self.global_xtx = np.zeros((n,m,m))
        for matrix in local_xt_matrices:
            local_xtx, local_xty = local_xt_matrices[matrix]
            self.global_xtx += local_xtx
            self.global_xty += local_xty
        return self.global_xtx, self.global_xty

