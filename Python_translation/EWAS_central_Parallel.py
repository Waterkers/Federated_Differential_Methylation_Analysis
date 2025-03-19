from pyexpat.errors import XML_ERROR_FEATURE_REQUIRES_XML_DTD
import numpy as np
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests
from tqdm.contrib.concurrent import process_map
import multiprocessing
import os
from functools import partial

def EWASCalculation(row, x_matrix):
    """
    Returns the regression coefficient, standard error, mu, sum of squares error, t statistic and p-value for
    each row (probe) of an y_matrix given the corresponding design matrix (x_matrix).

    """

    y_m = row
    x_t = x_matrix.T @ x_matrix
    x_t_y = x_matrix.T @ y_m
    x_t_inv = np.linalg.inv(x_t)
    coef = x_t_inv @ x_t_y  # beta
    #coefficient.append(coef)
    stan_er = np.sqrt(np.diag(x_t_inv))  # stdev_unscaled
    #standard_error.append(stan_er)
    m = x_matrix @ coef
    #mu.append(m)
    stersq = np.sum((y_m - m) ** 2)
    #SSE.append(stersq)

    t = coef / stan_er
    #t_stat.append(t)
    df = row.shape[0] - 2  # degrees of freedom is defined as number of observations - 2
    p = scipy.stats.t.sf(abs(t), df)
    #p_value.append(p)
    return coef, stan_er, m, stersq, t, p



def EWAS_central(design_matrix, beta_values):
    x_matrix = design_matrix.values
    y_matrix = beta_values.values


    n = y_matrix.shape[0] # select the number of rows of the beta matrix - #genes that the linear model will be calculated for
    m = x_matrix.shape[1] #select the number of columns from the design matrix

    corrected_pvalue = []

    # with multiprocessing.Pool(processes=max(os.cpu_count() - 2, 2)) as processPool:
    #     EWASCalculationPartial = partial(EWASCalculation, x_matrix=x_matrix)
    #     coefficient, standard_error, mu, SSE, t_stat, p_value = zip(*processPool.map(EWASCalculationPartial, y_matrix))
    EWASCalculationPartial = partial(EWASCalculation, x_matrix=x_matrix)
    coefficient, standard_error, mu, SSE, t_stat, p_value = zip(*process_map(EWASCalculationPartial, y_matrix))

    # for i in tqdm(range(0, n), token='7885595681:AAFmqfAkTqvAdOCRWF1RMl0s26XUmj5_Yv8', chat_id='5671488828'):
    #     y_m = y_matrix[i, :]
    #     x_t = x_matrix.T @ x_matrix
    #     x_t_y = x_matrix.T @ y_m
    #     x_t_inv = np.linalg.inv(x_t)
    #     coef = x_t_inv @ x_t_y #beta
    #     coefficient.append(coef)
    #     stan_er = np.sqrt(np.diag(x_t_inv)) #stdev_unscaled
    #     standard_error.append(stan_er)
    #     m = x_matrix @ coef
    #     mu.append(m)
    #     stersq = np.sum((y_m - mu)**2)
    #     SSE.append(stersq)
    #
    #     t = coef/stan_er
    #     t_stat.append(t)
    #     df = y_matrix.shape[1]-2 #degrees of freedom is defined as number of observations - 2
    #     p = scipy.stats.t.sf(abs(t), df)
    #     p_value.append(p)
    #
    #turn the results saved in lists into a dataframe for each covariate with the probe ids as index
    for i in range(0,len(p_value)):
        corrected_pvalue.append(multipletests(p_value[i], method="fdr_bh")[1])
    results_corp = pd.DataFrame(corrected_pvalue, index=beta_values.index, columns=design_matrix.columns)
    result_coef = pd.DataFrame(coefficient, index=beta_values.index, columns=design_matrix.columns)
    result_staner = pd.DataFrame(standard_error, index = beta_values.index, columns=design_matrix.columns)
    result_pvalue = pd.DataFrame(p_value, index=beta_values.index, columns=design_matrix.columns)

    #create a final results dataframe that contains the coefficient, standard error and p-value of the diagnosis covariate included in the linear regression
    results_EWAS = pd.concat([result_coef, result_staner, result_pvalue, results_corp], axis = 1, keys = ["Coefficient", "Standard Error", "P-value", "Corrected P-value"])
    return results_EWAS, SSE

if __name__ == '__main__':
    EWAS_central()