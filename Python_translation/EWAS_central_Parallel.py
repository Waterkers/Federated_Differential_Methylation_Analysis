import numpy as np
import pandas as pd
import scipy.stats
import yaml
from statsmodels.stats.multitest import multipletests
from tqdm.contrib.telegram import tqdm
import multiprocessing
import os
from functools import partial
import argparse


config_file = '/cosybio/project/vanElferen/FedEWAS/telegrambot_config.yml'
with open(config_file) as configStream:
    config = yaml.load(configStream, Loader=yaml.Loader)

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
    stan_er = np.sqrt(np.diag(x_t_inv))  # stdev_unscaled
    m = x_matrix @ coef
    stersq = np.sum((y_m - m) ** 2)
    t = coef / stan_er
    df = row.shape[0] - 2  # degrees of freedom is defined as number of observations - 2
    p = scipy.stats.t.sf(abs(t), df)
    return coef, stan_er, m, stersq, t, p#{'coef':coef, 'stan_er':stan_er, 'm':m, 'stersq':stersq, 't':t, 'p':p}



def EWAS_central(design_matrix, beta_values, identifier):
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
    with multiprocessing.Pool(processes=max(int(os.cpu_count()- 2), 2)) as processPool, tqdm(total=int(n),
                                                                                              token=
                                                                                              config[
                                                                                                  'serverProgressBot'][
                                                                                                  'token'],
                                                                                              chat_id=
                                                                                              config[
                                                                                                  'serverProgressBot'][
                                                                                                  'chat_id'], desc=f'Central EWAS {identifier}') as pbar:
        EWASCalculationPartial = partial(EWASCalculation, x_matrix=x_matrix)
        for result in processPool.imap(EWASCalculationPartial, y_matrix):
            # increment the progress bar
            pbar.update()
            pbar.refresh()
            print(len(result))
            coef, stan_er, m, stersq, t, p = result
            coefficient.append(coef)
            standard_error.append(stan_er)
            mu.append(m)
            SSE.append(stersq)
            t_stat.append(t)
            p_value.append(p)


    # EWASCalculationPartial = partial(EWASCalculation, x_matrix=x_matrix)
    # coefficient, standard_error, mu, SSE, t_stat, p_value = zip(*process_map(EWASCalculationPartial, y_matrix, max_workers=max(int(os.cpu_count() / 2), 2)))

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
    parser = argparse.ArgumentParser()
    parser.add_argument('designMatrix', help='design matrix')
    parser.add_argument('betaValues', help='beta values')
    EWAS_central(parser.parse_args().designMatrix, parser.parse_args().betaValues)
    # designMatrix = pd.read_csv('/home/silke/Documents/Fed_EWAS/Data/GSE134379_RAW/Small_EWAS_design.csv', index_col=0)
    # normalisedBetas = pd.read_csv(
    #     '/home/silke/Documents/Fed_EWAS/Data/QC_Python/GSE134379_half_normalised_betas_python.csv', index_col=0)
    # results_EWAS, SSE = EWAS_central(designMatrix, normalisedBetas.iloc[0:1000,:])
