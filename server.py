from threading import local
import sys
import pandas as pd
import numpy as np
import scipy
from statsmodels.stats.multitest import multipletests
#from scipy.stats import rankdata
#from scipy import interpolate
from scipy.special import digamma,polygamma
from scipy.stats import t
import re
import sys

def cov2cor(cov_coef):
    cor  = np.diag(cov_coef)**-0.5 * cov_coef
    cor = cor.T * np.diag(cov_coef)**-0.5
    np.fill_diagonal(cor,1)
    return cor

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
        
        self.cov_coef = None
        self.stdev_unscaled = None
        self.var = None
        self.sigma = None
        self.df_residual = None
        self.df_total = None
        self.results = None
        self.table = None
    
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
        global_conditions = self.variables.copy()
        global_conditions.extend(self.confounders)
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
            self.confounders.remove("Sentrix_ID")
        if self.global_PlateID:
            self.confounders.extend(self.global_PlateID)
            self.confounders.remove("Plate_ID")
        
        n = len(self.global_probes)
        m = (len(self.variables) + len(self.confounders))
        self.global_xtx = np.zeros((n,m,m))
        self.global_xty = np.zeros((n,m))
        self.stdev_unscaled = np.zeros((n,m))
        for i in range(0,len(local_xt_matrices)):
            local_xtx, local_xty = local_xt_matrices[i]
            self.global_xtx += local_xtx
            self.global_xty += local_xty

        self.beta = np.zeros((n,m))
        self.rank = np.ones(n)*m
    
        for i in range(0,n):
            invXtX  = np.linalg.inv(self.XtX_glob[i,:,:])
            self.beta[i,:] = invXtX @ self.XtY_glob[i,:]
            self.stdev_unscaled[i,:] = np.sqrt(np.diag(invXtX ))
        #return self.global_xtx, self.global_xty

    def aggregate_SSE_and_cov_coef(self,SSE_list,cov_coef_list, weighted=False):
        
        M= sum(self.samples_per_client) # total number of samples
        n = len(self.global_probes)
        k = len(self.confounders+self.variables)
        self.SSE = np.zeros(n)
        self.cov_coef = np.zeros((k,k))
        
        for c in range(0,len(self.client_names)):
            self.cov_coef += cov_coef_list[c]
            for i in range(0,n):
                self.SSE[i] += SSE_list[c][i]
            
        self.cov_coef = np.linalg.inv(self.cov_coef )        
        # estimated residual variance
        self.var = self.SSE/(M-k)
        # estimated residual standard deviations
        if weighted:
            self.sigma_w =  np.sqrt(self.var) 
        else:
            self.sigma =  np.sqrt(self.var)
        # degrees of freedom
        self.df_residual = np.ones(n)*(M-k)

    def make_contrasts(self, contrasts=[]):
        '''Creates contrast matrix given deisgn matrix and pairs or columns to compare.\n
        For example:\n
        contrasts = [([A],[B]),([A,B],[C,D])] defines two contrasts:\n
        A-B and (A and B) - (C and D).'''
        df = {}
        conditions = self.variables + self.confounders
        for contr in contrasts:
            group1 , group2 = contr
            for name in group1+group2:
                if not name in conditions:
                    print(name, "not found in the design matrix.",file=sys.stderr)
                    exit(1)
            contr_name = "".join(map(str,group1))+"_vs_"+"".join(map(str,group2))
            c=pd.Series(data=np.zeros(len(conditions)),index=conditions)
            c[group1] = 1
            c[group2] = -1
            df[contr_name] = c
        return (pd.DataFrame.from_dict(df))  

    
    def fit_contasts(self,contrast_matrix):
        ncoef = self.cov_coef.shape[1]
        #	Correlation matrix of estimable coefficients
        #	Test whether design was orthogonal
        if not np.any(self.cov_coef):
            print("no coef correlation matrix found in fit - assuming orthogonal",file=sys.stderr)
            cormatrix = np.identity(ncoef)
            orthog = True
        else:
            cormatrix = cov2cor(self.cov_coef)
            if cormatrix.shape[0]*cormatrix.shape[1] < 2: 
                orthog = True
            else:
                if np.sum(np.abs(np.tril(cormatrix,k=-1))) < 1e-12:
                    orthog = True
                else:
                    orthog = False
        #print("is design orthogonal:",orthog)
        #	Replace NA coefficients with large (but finite) standard deviations
        #	to allow zero contrast entries to clobber NA coefficients.
        if np.any(np.isnan(self.beta)):
            print("Replace NA coefficients with large (but finite) standard deviations",file=sys.stderr)
            np.nan_to_num(self.beta,nan=0)
            np.nan_to_num(self.stdev_unscaled, nan=1e30)

        self.beta = self.beta.dot(contrast_matrix)
        # New covariance coefficiets matrix
        self.cov_coef = contrast_matrix.T.dot(self.cov_coef).dot(contrast_matrix)

        if orthog:
            self.stdev_unscaled = np.sqrt((self.stdev_unscaled**2).dot(contrast_matrix**2))
        else:
            n_genes = self.beta.shape[0]
            U = np.ones((n_genes, contrast_matrix.shape[1])) # genes x contrasts
            o = np.ones(ncoef)
            R = np.linalg.cholesky(cormatrix).T
            for i in range(0,n_genes):
                RUC = R @ (self.stdev_unscaled[i,] * contrast_matrix.T).T
                U[i,] = np.sqrt(o @ RUC**2)
            self.stdev_unscaled = U

    #### e-Bayes ############

    def trigamma(self,x):
        return polygamma(1,x)

    def psigamma(self,x,deriv=2):
        return polygamma(deriv,x)

    def trigammaInverse(self,x):

        if not hasattr(x, '__iter__'):
            x_ = np.array([x])
        for i in range(0,x_.shape[0]):
            if np.isnan(x_[i]):
                x_[i]= np.NaN
            elif x>1e7:
                x_[i] = 1./np.sqrt(x[i])
            elif x< 1e-6:
                x_[i] = 1./x[i]
        # Newton's method
        y = 0.5+1.0/x_
        for i in range(0,50):
            tri = self.trigamma(y)
            dif = tri*(1.0-tri/x_)/self.psigamma(y,deriv=2)
            y = y+dif
            if(np.max(-dif/y) < 1e-8): # tolerance
                return y

        print("Warning: Iteration limit exceeded",file=sys.stderr)
        return y

    def fitFDist(self,x, df1, covariate=False):
        '''Given x (sigma^2) and df1 (df_residual), fits x ~ scale * F(df1,df2) and returns estimated df2 and scale (s0^2)'''
        if covariate:
            # TBD
            print("Set covariate=False.", file=sys.stderr)
            return
        # Avoid zero variances
        x = [max(x,0) for x in x]
        m = np.median(x)
        if(m==0):
            print("Warning: More than half of residual variances are exactly zero: eBayes unreliable", file=sys.stderr)
            m = 1
        else:
            if 0 in x: 
                print("Warning: Zero sample variances detected, have been offset (+1e-5) away from zero", file=sys.stderr)

        x = [max(x,1e-5 * m) for x in x] 
        z = np.log(x)
        e = z-digamma(df1*1.0/2)+np.log(df1*1.0/2)
        emean = np.nanmean(e)
        evar = np.nansum((e-emean)**2)/(len(x)-1)

        # Estimate scale and df2
        evar = evar - np.nanmean(self.trigamma(df1*1.0/2))
        
        if evar > 0:
            df2 = 2*self.trigammaInverse(evar)
            s20 = np.exp(emean+digamma(df2*1.0/2)-np.log(df2*1.0/2))
        else:
            df2 = np.Inf
            s20 = np.exp(emean)

        return s20,df2

    def posterior_var(self, var_prior=np.ndarray([]), df_prior=np.ndarray([])):
        var = self.var
        df=self.df_residual
        ndxs = np.argwhere(np.isfinite(var)).reshape(-1)
        # if not infinit vars
        if len(ndxs)==len(var): #np.isinf(df_prior).any():
            return (df*var + df_prior*var_prior) / (df+df_prior) # var_post  
        #For infinite df.prior, set var_post = var_prior
        var_post = np.repeat(var_prior,len(var))
        for ndx in ndxs:
            var_post[ndx] = (df[ndx]*var[ndx] + df_prior*var_prior)/(df[ndx]+df_prior)
        return var_post

    def squeezeVar(self, covariate=False, robust=False, winsor_tail_p=(0.05,0.1)):
        '''Estimates df and var priors and computes posterior variances.'''
        if robust:
            # TBD fitFDistRobustly()
            print("Set robust=False.",file=sys.stderr)
            return
        else:
            var_prior, df_prior = self.fitFDist(self.var, self.df_residual, covariate=covariate)

        if np.isnan(df_prior):
            print ("Error: Could not estimate prior df.",file=sys.stderr)
            return

        var_post = self.posterior_var(var_prior=var_prior, df_prior=df_prior)
        self.results = {"df_prior":df_prior,"var_prior":var_prior,"var_post":var_post}

    def moderatedT(self,covariate=False,robust=False, winsor_tail_p=(0.05,0.1)):
        #var,df_residual,coefficients,stdev_unscaled,
        self.squeezeVar(covariate=covariate, robust=robust, winsor_tail_p=winsor_tail_p)
        
        self.results["s2_prior"] = self.results["var_prior"]
        self.results["s2_post"] = self.results["var_post"]
        del self.results["var_prior"]
        del self.results["var_post"]
        self.results["t"] = self.beta / self.stdev_unscaled
        self.results["t"] = self.results["t"].T / np.sqrt( self.results["s2_post"])
        self.df_total = self.df_residual + self.results["df_prior"]
        df_pooled = sum(self.df_residual)
        self.df_total = np.minimum(self.df_total,df_pooled) # component-wise min

        self.results["p_value"] = 2*t.cdf(-np.abs(self.results["t"]),df=self.df_total)
        self.results["p_value"] = self.results["p_value"].T
        self.results["t"] = self.results["t"].T
        return self.results

    def tmixture_matrix(self,var_prior_lim=False,proportion=0.01):
        tstat = self.results["t"]
        stdev_unscaled = self.stdev_unscaled
        df_total = self.df_total
        ncoef = self.results["t"].shape[1]
        v0 = np.zeros(ncoef)
        for j in range(0,ncoef):
            v0[j] = self.tmixture_vector(tstat[:,j],stdev_unscaled[:,j],df_total,proportion,var_prior_lim)
        return v0

    def tmixture_vector(self,tstat,stdev_unscaled,df,proportion,var_prior_lim):
        ngenes = len(tstat)

        #Remove missing values
        notnan_ndx = np.where(~np.isnan(tstat))[0]
        if len(notnan_ndx) < ngenes:
            tstat = tstat[notnan_ndx]
            stdev_unscaled = stdev_unscaled[notnan_ndx]
            df = df[notnan_ndx]

        # ntarget t-statistics will be used for estimation

        ntarget = int(np.ceil(proportion/2*ngenes))
        if ntarget < 1: #
            return 

        # If ntarget is v small, ensure p at least matches selected proportion
        # This ensures ptarget < 1
        p = np.maximum(ntarget*1.0/ngenes,proportion)

        #Method requires that df be equal
        tstat = abs(tstat)
        MaxDF = np.max(df)
        i = np.where( df < MaxDF)[0]
        if len(i)>0:
            TailP = t.logcdf(tstat[i],df=df[i]) # PT: CDF of t-distribution: pt(tstat[i],df=df[i],lower.tail=FALSE,log.p=TRUE)
            # QT - qunatile funciton - returns a threshold value x 
            # below which random draws from the given CDF would fall p percent of the time. [wiki]
            tstat[i] = t.ppf(np.exp(TailP),df=MaxDF) # QT: qt(TailP,df=MaxDF,lower.tail=FALSE,log.p=TRUE)
            df[i] = MaxDF

        #Select top statistics
        order = tstat.argsort()[::-1][:ntarget] # TBD: ensure the order is decreasing
        tstat = tstat[order]
        v1 =  stdev_unscaled[order]**2


        #Compare to order statistics
        rank = np.array(range(1,ntarget+1))
        p0 =  2*t.sf(tstat,df=MaxDF) # PT
        ptarget = ((rank-0.5)/ngenes - (1.0-p)*p0) / p
        v0 = np.zeros(ntarget)
        pos = np.where(ptarget > p0)[0]
        if len(pos)>0:
            qtarget = -t.ppf(ptarget[pos]/2,df=MaxDF) #qt(ptarget[pos]/2,df=MaxDF,lower.tail=FALSE)
            #print(qtarget[:5])
            v0[pos] = v1[pos]*((tstat[pos]/qtarget)**2-1)

        if var_prior_lim[0] and var_prior_lim[1]:
            v0 = np.minimum(np.maximum(v0,var_prior_lim[0]),var_prior_lim[1])

        return np.mean(v0)
    
    
    def Bstat(self,stdev_coef_lim = np.array([0.1,4]),proportion = 0.01):

        var_prior_lim  = stdev_coef_lim**2/np.median(self.results["s2_prior"])
        #print("Limits for var.prior:",var_prior_lim)

        self.results["var_prior"] = self.tmixture_matrix(proportion=0.01,var_prior_lim=var_prior_lim)

        nan_ndx = np.argwhere(np.isnan(self.results["var_prior"]))
        if len(nan_ndx)>0:
            self.results["var.prior"][ nan_ndx] <- 1.0/self.results["s2_prior"]
            print("Warning: Estimation of var.prior failed - set to default value",file = sys.stderr)
        r = np.outer(np.ones(self.results["t"].shape[0]),self.results["var_prior"])
        r = (self.stdev_unscaled**2+r) / self.stdev_unscaled**2
        t2 = self.results["t"]**2

        valid_df_ndx = np.where(self.results["df_prior"] <= 1e6)[0]
        if len(valid_df_ndx)<len(self.results["df_prior"]):
            print("Large (>1e6) priors for DF:", len(valid_df_ndx))
            kernel = t2*(1-1.0/r)/2
            for i in valid_df_ndx:
                kernel[i] = (1+self.df_total[i])/2*np.log((t2[i,:].T+self.df_total[i]) / ((t2[i,:]/r[i,:]).T+self.df_total[i]))
        else:
            kernel = (1+self.df_total)/2*np.log((t2.T+self.df_total)/((t2/r).T+self.df_total))

        self.results["lods"] = np.log(proportion/(1-proportion))-np.log(r)/2+kernel.T
    
    def topTableT(self,adjust="fdr_bh", p_value=1.0,lfc=0,confint=0.95):
        feature_names = self.global_genes
        self.results["logFC"] = pd.Series(self.beta[:,0],index=feature_names)
        
        # confidence intervals for LogFC
        if confint:
            alpha = (1.0+confint)/2
            margin_error = np.sqrt(self.results["s2_post"]) *self.stdev_unscaled[:,0] * t.ppf(alpha, df=self.df_total)
            self.results["CI.L"] = self.results["logFC"]-margin_error
            self.results["CI.R"] = self.results["logFC"] +margin_error
        # adjusting p-value for multiple testing
        if_passed, adj_pval,alphacSidak,alphacBonf = multipletests(self.results["p_value"][:,0], alpha=p_value, method=adjust,
                                           is_sorted=False, returnsorted=False)
        self.results["adj.P.Val"] = pd.Series(adj_pval,index=feature_names)
        self.results["P.Value"] = pd.Series(self.results["p_value"][:,0],index=feature_names)
        # make table 
        self.table = np.copy(self.results)
        # remove 'df_prior', 's2_prior', 's2_post', 'df_total','var_prior'
        for key in ['df_prior', 's2_prior', 's2_post', 'var_prior',"p_value"]:        
            del self.table[key]
        self.table["t"] = pd.Series(self.table["t"][:,0],index=feature_names)
        self.table["lods"] = pd.Series(self.table["lods"][:,0],index=feature_names)
        self.table = pd.DataFrame.from_dict(self.table)
    
    def eBayes(self):
        covariate = False # Amean for limma-trend
        robust = False # 
        winsor_tail_p = (0.05,0.1) # needed for fitFDistRobustly()

        var = self.sigma**2 
        self.results = self.moderatedT(covariate=covariate,robust=robust, winsor_tail_p=winsor_tail_p)
        #self.results = moderatedT(self.var,self.df_residual,self.beta,self.stdev_unscaled,
        #                         covariate=covariate,robust=robust, winsor_tail_p=winsor_tail_p)
        self.results["AveExpr"] = self.Amean
        
        self.Bstat(stdev_coef_lim = np.array([0.1,4]),proportion = 0.01)
        
        self.topTableT(adjust="fdr_bh", p_value=1.0,lfc=0,confint=0.95)
        self.table = self.table.sort_values(by="P.Value")
        


""" 
    def calculate_EWAS_results(self):
        n = len(self.global_probes)
        m = (len(self.variables) + len(self.confounders))
        coef = np.zeros((n,m))
        stnd_err = np.zeros((n,m))
        p_value = np.zeros((n,m))
        p_value_cor = np.zeros((n,m))
        mchange = np.zeros((n,m))
        
        for i in range(0,n):
            xtx_inv = np.linalg.inv(self.global_xtx[i,:,:])
            coef[i,:] = xtx_inv @ self.global_xty[i,:]
            stnd_err[i,:] = np.sqrt(np.diag(xtx_inv))
            t = coef[i,:]/stnd_err[i,:]
            df = m-2
            #for j in range(0,m):
            p_value[i,:] = scipy.stats.t.sf(t, df)#*2 # using the absolute t value
            #self.p_value_cor[i,:] = multipletests(self.p_value[i,:], method="fdr_bh")[1]
        # create EWAS results dataframe with all information grouped by variable/confounder
        for i in range(0,m):
            p_value_cor[:,i] = multipletests(p_value[:,i], method="fdr_bh")[1]
        
        for i in range(0,m):
            mchange[:,i] = coef[:,i]
        #coef = pd.DataFrame(self.coef,index=self.probes, columns= self.designcolumns)
        #stnErr = pd.DataFrame(self.stnd_err, index=self.probes, columns= self.designcolumns)
        #p_val = pd.DataFrame(self.p_value, index=self.probes, columns= self.designcolumns)
        #p_val_corrected = pd.DataFrame(self.p_value_cor, index=self.probes, columns= self.designcolumns)
        #mchange = pd.DataFrame(self.mchange, index = self.probes, columns=self.designcolumns)
        # create a dataframe with the corrected p-values
        
        self.EWAS_results = pd.DataFrame({"Coefficient":coef[:,0], "Standard Error":stnd_err[:,0], 
            "P-value":p_value[:,0], "Corrected P-value":p_value_cor[:,0], "Methylation Change": mchange[:,0]},
            index=self.global_probes)
        
        return self.EWAS_results """

