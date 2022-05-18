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
        self.cohort_effects = []
        self.global_SentrixID = []
        self.global_PlateID = []

        self.global_xtx = None
        self.global_xty = None
        
    
    def get_clients(self, cohort_name, client_probes, client_n_samples):
        '''
        get client information
        '''
        self.client_names.append(cohort_name)
        self.samples_per_client.append(len(client_n_samples))
        self.client_total_probes.append(client_probes)
        self.cohort_effects = sorted(self.client_names)[:-1]
        
    def find_cohort_effects(self):
        self.cohort_effects = sorted(self.client_names)[:-1]
        return self.cohort_effects
           
    def find_global_probes(self):
        '''
        find the probes that are present in the data held by the clients 
        '''
        global_probes = set(self.client_total_probes[0])
        for probe_list in self.client_total_probes[1:]:
            global_probes = global_probes.intersection(set(probe_list))
        self.global_probes = list(global_probes)
        return self.global_probes
        """ 
        n_clients = len(self.client_names)
        #common_index = []
        for probe in self.client_total_probes:
            if (self.client_total_probes.count(probe) == n_clients):
                self.global_probes.append(probe)
            else: 
                pass  
        return self.global_probes
 """
    def return_global_conditions(self):
        '''
        send the global conditions (variables and confounders) included in the analysis to the client
        '''
        self.confounders.extend(self.variables)
        global_conditions = self.confounders
        return global_conditions

    def return_global_SentrixID(self, *local_SentrixID):
        '''
        send all but one sentrixIDs present across all clients back to each client so they can be included as
        confounders in the design matrix

        This function will only be used if SentrixID is inlcuded in the global conditions when initialising
        the project.  
        '''
        sentrix_ID = local_SentrixID[0]
        for ids in local_SentrixID[1:]:
            sentrix_ID.extend(ids)
        self.global_SentrixID = list(set(sentrix_ID))[:-1]
        
        return self.global_SentrixID

    def return_global_PlateID(self, *local_PlateID):
        '''
        send all but one plateIDs present across all clients back to each client so they can be included as
        confounders in the design matrix

        This function will only be used if PlateID is inlcuded in the global conditions when initialising
        the project.  
        '''
        plate_ID = local_PlateID[0]
        for ids in local_PlateID[1:]:
            plate_ID.extend(ids)
        self.global_PlateID = list(set(plate_ID))[:-1]
        
        return self.global_PlateID
           
    
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
    
    
    def global_regression_parameter(self, *local_xt_matrices):
        '''
        aggregate local xt_x to global xt_x matrix for linear regression
        '''
        
        self.confounders.extend(self.cohort_effects)
        if self.global_SentrixID:
            self.confounders.extend(self.global_SentrixID)
        if self.global_PlateID:
            self.confounders.extend(self.global_PlateID)
        n = len(self.global_probes)
        m = (len(self.variables) + len(self.confounders))
        self.global_xtx = np.zeros((n,m,m))
        self.global_xty = np.zeros((n,m))
        for i in range(0,len(local_xt_matrices)):
            local_xtx, local_xty = local_xt_matrices[i]
            self.global_xtx += local_xtx
            self.global_xty += local_xty
        return self.global_xtx, self.global_xty

