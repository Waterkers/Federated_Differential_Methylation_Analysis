from pyexpat.errors import XML_ERROR_FEATURE_REQUIRES_XML_DTD


def EWAS_central(design_matrix, beta_values):
    import numpy as np
    import pandas as pd
    import scipy.stats
    from statsmodels.stats.multitest import multipletests
    x_matrix = design_matrix.values
    y_matrix = beta_values.values


    n = y_matrix.shape[0] # select the number of rows of the beta matrix - #genes that the linear model will be calculated for
    m = x_matrix.shape[1] #select the number of columns from the design matrix

    

    coefficient = []
    standard_error = []
    mu = []
    SSE = []
    t_stat = []
    p_value = []
    corrected_pvalue = []

    for i in range(0, n):
        y_m = y_matrix[i, :]
        x_t = x_matrix.T @ x_matrix
        x_t_y = x_matrix.T @ y_m
        x_t_inv = np.linalg.inv(x_t)
        coef = x_t_inv @ x_t_y #beta
        coefficient.append(coef)
        stan_er = np.sqrt(np.diag(x_t_inv)) #stdev_unscaled
        standard_error.append(stan_er)
        m = x_matrix @ coef
        mu.append(m)
        stersq = np.sum((y_m - mu)**2)
        SSE.append(stersq)
        
        t = coef/stan_er
        t_stat.append(t)
        df = y_matrix.shape[1]-2 #degrees of freedom is defined as number of observations - 2 
        p = scipy.stats.t.sf(abs(t), df)
        p_value.append(p)
    
    #turn the results saved in lists into a dataframe for each covariate with the probe ids as index
    for i in range(0,len(p_value)):
        corrected_pvalue.append(multipletests(p_value[i], method="fdr_bh")[1])
    results_corp = pd.DataFrame(corrected_pvalue, index=beta_values.index, columns=design_matrix.columns)
    result_coef = pd.DataFrame(coefficient, index=beta_values.index, columns=design_matrix.columns)
    result_staner = pd.DataFrame(standard_error, index = beta_values.index, columns=design_matrix.columns)
    result_pvalue = pd.DataFrame(p_value, index=beta_values.index, columns=design_matrix.columns)

    #create a final results dataframe that contains the coefficient, standard error and p-value of the diagnosis covariate included in the linear regression
    results_EWAS = pd.concat([result_coef, result_staner, result_pvalue, results_corp], axis = 1, keys = ["Coefficient", "Standard Error", "P-value", "Corrected P-value"])
    return results_EWAS