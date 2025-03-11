import numpy as np
import pandas as pd
from copy import copy
from scipy.special import digamma,polygamma
import sys
from scipy.stats import t
from statsmodels.stats.multitest import multipletests


class eBayesLocal:
    def __init__(self,df_path, sse_path, n_samples):
        df = df_path # pd.read_csv(df_path, index_col=0, header=[0,1])
        #print(df.head())
        self.linearModel = df
        self.n_samples = n_samples
        self.beta = df.loc[:, pd.IndexSlice['Coefficient',:]].droplevel(0, axis=1)
        self.stdev_unscaled = df.loc[:, pd.IndexSlice['Standard Error',:]].droplevel(0, axis=1)
        self.SSE = sse_path
        #with open(sse_path, 'r') as file:
        #    for i in file.readlines():
        #        self.SSE.append(i)
        self.SSE = np.array(self.SSE)
        if (self.n_samples-self.beta.shape[1]) == 0:
            print('degrees of freedom technically 0 because too few samples, for debug purposes set degrees of freedom to 1')
            self.df_residual = np.ones(self.beta.shape[0]) * 1
            self.var = self.SSE / 1
        else:
            self.df_residual = np.ones(self.beta.shape[0]) * (self.n_samples-self.beta.shape[1])
            self.var = self.SSE / (self.n_samples - self.beta.shape[1])
        self.table = None
        self.sigma = np.sqrt(self.var)

    def trigamma(self ,x):
        return polygamma(1 ,x)

    def psigamma(self ,x ,deriv=2):
        return polygamma(deriv ,x)

    def trigammaInverse(self ,x):

        if not hasattr(x, '__iter__'):
            x_ = np.array([x])
        for i in range(0 ,x_.shape[0]):
            if np.isnan(x_[i]):
                x_[i ]= np.NaN
            elif x> 1e7:
                x_[i] = 1. / np.sqrt(x[i])
            elif x < 1e-6:
                x_[i] = 1. / x[i]
        # Newton's method
        y = 0.5 + 1.0 / x_
        for i in range(0, 50):
            tri = self.trigamma(y)
            dif = tri * (1.0 - tri / x_) / self.psigamma(y, deriv=2)
            y = y + dif
            if (np.max(-dif / y) < 1e-8):  # tolerance
                return y

        print("Warning: Iteration limit exceeded", file=sys.stderr)
        return y


    def fitFDist(self, x, df1, covariate=False):
        '''Given x (sigma^2) and df1 (df_residual), fits x ~ scale * F(df1,df2) and returns estimated df2 and scale (s0^2)'''
        if covariate:
            # TBD
            print("Set covariate=False.", file=sys.stderr)
            return
        # Avoid zero variances
        x = [max(x, 0) for x in x]
        m = np.median(x)
        if (m == 0):
            print("Warning: More than half of residual variances are exactly zero: eBayes unreliable", file=sys.stderr)
            m = 1
        else:
            if 0 in x:
                print("Warning: Zero sample variances detected, have been offset (+1e-5) away from zero", file=sys.stderr)

        x = [max(x, 1e-5 * m) for x in x]
        z = np.log(x)
        e = z - digamma(df1 * 1.0 / 2) + np.log(df1 * 1.0 / 2)
        emean = np.nanmean(e)
        evar = np.nansum((e - emean) ** 2) / (len(x) - 1)

        # Estimate scale and df2
        evar = evar - np.nanmean(self.trigamma(df1 * 1.0 / 2))

        if evar > 0:
            df2 = 2 * self.trigammaInverse(evar)
            s20 = np.exp(emean + digamma(df2 * 1.0 / 2) - np.log(df2 * 1.0 / 2))
        else:
            df2 = np.Inf
            s20 = np.exp(emean)

        return s20, df2


    def posterior_var(self, var_prior=np.ndarray([]), df_prior=np.ndarray([])):
        var = self.var
        df = self.df_residual
        ndxs = np.argwhere(np.isfinite(var)).reshape(-1)
        # if not infinit vars
        if len(ndxs) == len(var):  # np.isinf(df_prior).any():
            return (df * var + df_prior * var_prior) / (df + df_prior)  # var_post
        # For infinite df.prior, set var_post = var_prior
        var_post = np.repeat(var_prior, len(var))
        for ndx in ndxs:
            var_post[ndx] = (df[ndx] * var[ndx] + df_prior * var_prior) / (df[ndx] + df_prior)
        return var_post


    def squeezeVar(self, covariate=False, robust=False, winsor_tail_p=(0.05, 0.1)):
        '''Estimates df and var priors and computes posterior variances.'''
        if robust:
            # TBD fitFDistRobustly()
            print("Set robust=False.", file=sys.stderr)
            return
        else:
            var_prior, df_prior = self.fitFDist(self.var, self.df_residual, covariate=covariate)

        if np.isnan(df_prior):
            print("Error: Could not estimate prior df.", file=sys.stderr)
            return

        var_post = self.posterior_var(var_prior=var_prior, df_prior=df_prior)
        self.results = {"df_prior": df_prior, "var_prior": var_prior, "var_post": var_post}


    def moderatedT(self, covariate=False, robust=False, winsor_tail_p=(0.05, 0.1)):
        # var,df_residual,coefficients,stdev_unscaled,
        self.squeezeVar(covariate=covariate, robust=robust, winsor_tail_p=winsor_tail_p)

        self.results["s2_prior"] = self.results["var_prior"]
        self.results["s2_post"] = self.results["var_post"]
        del self.results["var_prior"]
        del self.results["var_post"]
        self.results["t"] = self.beta / self.stdev_unscaled
        self.results["t"] = self.results["t"].T / np.sqrt(self.results["s2_post"])
        self.df_total = self.df_residual + self.results["df_prior"]
        df_pooled = sum(self.df_residual)
        self.df_total = np.minimum(self.df_total, df_pooled)  # component-wise min

        self.results["p_value"] = 2 * t.cdf(-np.abs(self.results["t"]), df=self.df_total)
        self.results["p_value"] = self.results["p_value"].T
        self.results["t"] = self.results["t"].T
        return self.results


    def tmixture_matrix(self, var_prior_lim=False, proportion=0.01):
        tstat = self.results["t"]
        stdev_unscaled = self.stdev_unscaled
        df_total = self.df_total
        ncoef = self.results["t"].shape[1]
        v0 = np.zeros(ncoef)
        for j in range(0, ncoef):
            v0[j] = self.tmixture_vector(tstat[:, j], stdev_unscaled[:, j], df_total, proportion, var_prior_lim)
        return v0


    def tmixture_vector(self, tstat, stdev_unscaled, df, proportion, var_prior_lim):
        ngenes = len(tstat)

        # Remove missing values
        notnan_ndx = np.where(~np.isnan(tstat))[0]
        if len(notnan_ndx) < ngenes:
            tstat = tstat[notnan_ndx]
            stdev_unscaled = stdev_unscaled[notnan_ndx]
            df = df[notnan_ndx]

        # ntarget t-statistics will be used for estimation

        ntarget = int(np.ceil(proportion / 2 * ngenes))
        if ntarget < 1:  #
            return

            # If ntarget is v small, ensure p at least matches selected proportion
        # This ensures ptarget < 1
        p = np.maximum(ntarget * 1.0 / ngenes, proportion)

        # Method requires that df be equal
        tstat = abs(tstat)
        MaxDF = np.max(df)
        i = np.where(df < MaxDF)[0]
        if len(i) > 0:
            TailP = t.logcdf(tstat[i],
                             df=df[i])  # PT: CDF of t-distribution: pt(tstat[i],df=df[i],lower.tail=FALSE,log.p=TRUE)
            # QT - qunatile funciton - returns a threshold value x
            # below which random draws from the given CDF would fall p percent of the time. [wiki]
            tstat[i] = t.ppf(np.exp(TailP), df=MaxDF)  # QT: qt(TailP,df=MaxDF,lower.tail=FALSE,log.p=TRUE)
            df[i] = MaxDF

        # Select top statistics
        order = tstat.argsort()[::-1][:ntarget]  # TBD: ensure the order is decreasing
        tstat = tstat[order]
        v1 = stdev_unscaled[order] ** 2

        # Compare to order statistics
        rank = np.array(range(1, ntarget + 1))
        p0 = 2 * t.sf(tstat, df=MaxDF)  # PT
        ptarget = ((rank - 0.5) / ngenes - (1.0 - p) * p0) / p
        v0 = np.zeros(ntarget)
        pos = np.where(ptarget > p0)[0]
        if len(pos) > 0:
            qtarget = -t.ppf(ptarget[pos] / 2, df=MaxDF)  # qt(ptarget[pos]/2,df=MaxDF,lower.tail=FALSE)
            # print(qtarget[:5])
            v0[pos] = v1[pos] * ((tstat[pos] / qtarget) ** 2 - 1)

        if var_prior_lim[0] and var_prior_lim[1]:
            v0 = np.minimum(np.maximum(v0, var_prior_lim[0]), var_prior_lim[1])

        return np.mean(v0)


    def Bstat(self, stdev_coef_lim=np.array([0.1, 4]), proportion=0.01):
        var_prior_lim = stdev_coef_lim ** 2 / np.median(self.results["s2_prior"])
        # print("Limits for var.prior:",var_prior_lim)

        self.results["var_prior"] = self.tmixture_matrix(proportion=0.01, var_prior_lim=var_prior_lim)

        nan_ndx = np.argwhere(np.isnan(self.results["var_prior"]))
        if len(nan_ndx) > 0:
            self.results["var.prior"][nan_ndx] < - 1.0 / self.results["s2_prior"]
            print("Warning: Estimation of var.prior failed - set to default value", file=sys.stderr)
        r = np.outer(np.ones(self.results["t"].shape[0]), self.results["var_prior"])
        r = (self.stdev_unscaled ** 2 + r) / self.stdev_unscaled ** 2
        t2 = self.results["t"] ** 2

        valid_df_ndx = np.where(self.results["df_prior"] <= 1e6)[0]
        if len(valid_df_ndx) < len(self.results["df_prior"]):
            print("Large (>1e6) priors for DF:", len(valid_df_ndx))
            kernel = t2 * (1 - 1.0 / r) / 2
            for i in valid_df_ndx:
                kernel[i] = (1 + self.df_total[i]) / 2 * np.log(
                    (t2[i, :].T + self.df_total[i]) / ((t2[i, :] / r[i, :]).T + self.df_total[i]))
        else:
            kernel = (1 + self.df_total) / 2 * np.log((t2.T + self.df_total) / ((t2 / r).T + self.df_total))

        self.results["lods"] = np.log(proportion / (1 - proportion)) - np.log(r) / 2 + kernel.T


    def topTableT(self, adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95):
        feature_names = self.beta.index.to_list()
        self.results["logFC"] = pd.Series(self.beta[:, 0], index=feature_names)

        # confidence intervals for LogFC
        if confint:
            alpha = (1.0 + confint) / 2
            margin_error = np.sqrt(self.results["s2_post"]) * self.stdev_unscaled[:, 0] * t.ppf(alpha, df=self.df_total)
            self.results["CI.L"] = self.results["logFC"] - margin_error
            self.results["CI.R"] = self.results["logFC"] + margin_error
        # adjusting p-value for multiple testing
        if_passed, adj_pval, alphacSidak, alphacBonf = multipletests(self.results["p_value"][:, 0], alpha=p_value,
                                                                     method=adjust,
                                                                     is_sorted=False, returnsorted=False)
        self.results["adj.P.Val"] = pd.Series(adj_pval, index=feature_names)
        self.results["P.Value"] = pd.Series(self.results["p_value"][:, 0], index=feature_names)
        # make table
        self.table = copy(self.results)
        # remove 'df_prior', 's2_prior', 's2_post', 'df_total','var_prior'
        for key in ['df_prior', 's2_prior', 's2_post', 'var_prior', "p_value"]:
            del self.table[key]
        self.table["t"] = pd.Series(self.table["t"][:, 0], index=feature_names)
        self.table["lods"] = pd.Series(self.table["lods"][:, 0], index=feature_names)
        self.table = pd.DataFrame.from_dict(self.table)


    def eBayes(self):
        covariate = False  # Amean for limma-trend
        robust = False  #
        winsor_tail_p = (0.05, 0.1)  # needed for fitFDistRobustly()

        var = self.sigma ** 2
        self.results = self.moderatedT(covariate=covariate, robust=robust, winsor_tail_p=winsor_tail_p)
        # self.results = moderatedT(self.var,self.df_residual,self.beta,self.stdev_unscaled,
        #                         covariate=covariate,robust=robust, winsor_tail_p=winsor_tail_p)
        # self.results["AveExpr"] = self.Amean

        self.Bstat(stdev_coef_lim=np.array([0.1, 4]), proportion=0.01)

        self.topTableT(adjust="fdr_bh", p_value=1.0, lfc=0, confint=0.95)
        self.table = self.table.sort_values(by="P.Value")