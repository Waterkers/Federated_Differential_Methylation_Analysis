def clusterGender(betas, sex, npcs = 20, thres = 0.5, makePlot = True, inLine = False):
    
    from sklearn.decomposition import PCA 
    import numpy as np
    import pandas as pd
    from scipy.stats import pearsonr
    from sklearn.cluster import KMeans
    from matplotlib import pyplot as plt

    random_state = 42 #to make sure to PCA output is stabel when the code is run multiple times on the same dataset

    if npcs > len(betas.columns):
        npcs = len(betas.columns)
        print("As there are only", len(betas.columns), "samples, this is the total number of PCs that can be looked at")
    
    
    complete_betas = betas.dropna(axis = 0, how = "any") # drop any rows - probes - that contain missing/NaN values
    pca = PCA(n_components = npcs)
    pca_results = pca.fit(complete_betas) #also performs zero-centering before fitting the PCA model, same as the r function used in the original code
    
    pca_iterator = np.arange(0, npcs)
    pca_cor = []
    codes_sex, options_sex = pd.factorize(sex)
    for id, x in np.ndenumerate(pca_iterator):
        cor, p_value = pearsonr(pca_results.components_[x,], codes_sex)
        pca_cor.append(float(cor)) # turn into list append
    

    pc_output = pd.DataFrame({"PC": pca_iterator, "Correlation": pca_cor})
    pc_output.sort_values(by = "Correlation", key = pd.Series.abs, ascending = False, inplace=True) #sort the correlation values in descending order based on their absolute value
    first = abs(pc_output.iloc[0]) #get the two highest correlated PCs
    second = abs(pc_output.iloc[1])
    print("The top correlated principle components with sex:", first.name + 1, second.name + 1) #the +1 is necessary to get the correct PC number when compared to R since python
    # indexes from 0 and R from 1
    
    """ if makePlot == True:
        # insert code for correlation plot
        if inLine == True:
            %matplotlib inline
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])
            ax.scatter(pca_results.components_[np.int_(first[0])], pca_results.components_[np.int_(second[0])], c = codes_sex, cmap="Set2")
            
        else
            fig = plt.figure()
            ax = fig.add_axes([0,0,1,1])
            ax.scatter(pca_results.components_[np.int_(first[0])], pca_results.components_[np.int_(second[0])], c = codes_sex, cmap="Set2")
      """ #needs to be fixed so the legend works properly  
        
    
    # Lousis kmean clustering algorithm with hartigan approximation I think
    KMeans_sex = KMeans(n_clusters = 2, random_state = random_state).fit(pca_results.components_[pc_output[abs(pc_output) > 0.2].dropna().index,].transpose())
    predSex_num = KMeans_sex.labels_
    
    #generate a dataframe or numpy array (whichever turns out to be easier) that contains the preditect sex based on the kmeans clustering labels    
    predSex = np.asarray(predSex_num).astype("object")
       
    maxSample = np.argmax(abs(pca_results.components_[first.name]))
    SamplePredSex = predSex_num[maxSample]
    SampleRealSex = sex[maxSample]
    SampleSex = np.asarray(options_sex == SampleRealSex).nonzero()

    Sex_sample = np.str(options_sex[SampleSex[0]].tolist()[0])
    Sex_other = np.str(options_sex[~SampleSex[0]].tolist()[0])

    predSex = np.where(predSex == SamplePredSex, Sex_sample, Sex_other)
    
    return predSex
    
    #hartigans kmeans clustering (fortran knms integration) in case the above yields deviating results
    