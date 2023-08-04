#functions used in TRAP analysis

import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt

def ssmd(df1, df2, sf=0.01):
    """This function calculates the ssmd
    https://en.wikipedia.org/wiki/Strictly_standardized_mean_difference
    Variables:
        df1 : Pandas dataframe
        df2 : Pandas dataframe 
        sf : smoothing factor"""
    ssmd = (df1.mean(axis=1) - df2.mean(axis=1)) / (( sf + (df1.std(axis=1)**2) + (df2.std(axis=1)**2) )**0.5)
    return ssmd



def log2_fold_change(df1, df2, averaging='mean', sf=0.01):
    """This function calculates the fold change
    Variables:
        df1 : Pandas dataframe
        df2 : Pandas dataframe
        averaging : str (mean or median)
        """
    if averaging == 'mean':
        fc = (df1.mean(axis=1)+sf)/(df2.mean(axis=1)+sf)
    elif averaging == 'median':
        fc = (df1.median(axis=1)+sf)/(df2.median(axis=1)+sf)

    fc= np.log2(fc)
    return fc


#functions used in TRAP analysis



def bhatt_distance(df1, df2, sf=0.01, sign_adjust=True):
    """Replicated from https://en.wikipedia.org/wiki/Bhattacharyya_distance
    where df1 : pandas data frame, experimental group  (p)
        df2 : pandas data frame, control group (q)
        sf : smoothing factor"""

    df1std = df1.std(axis=1)
    df2std = df2.std(axis=1)
    df1mean = df1.mean(axis=1)
    df2mean = df2.mean(axis=1)
    
    #Calculate a1=vardf1/ vardf2, a2=vardf2/vardf1
    a1 = (df1std)**2 / ((df2std**2).clip(lower=sf))
    a2 = (df2std)**2 / ((df1std**2).clip(lower=sf))
    a = 0.25*np.log(0.5+0.25*(a1 + a2)) #if we don't specify base, log is ln in numpy
    
    
    #b is the second part, with b1 as the numerator and b2 as the denominator
    b1 = (df1mean - df2mean)**2
    # b2 = (df1std+sf)**2 + (df2std+sf)**2
    b2 = (df1std**2 + df2std**2).clip(lower=sf)

    b = 0.25*(b1 / b2)

    #Summate to calc DB
    DB = a+b

    #calculate BC from DB
    BC = (1- (np.exp(-DB)).clip(upper=1) )

    if sign_adjust is True:
        fc_sign = np.sign(log2_fold_change(df1, df2))
        #Multiply by sign
        BC = BC*fc_sign


    return BC


def kdensity(df1, sf=0.01):

    fig,ax=plt.subplots(figsize=(8, 6))
    for i in range(0,len(df1.columns)):
        P1 = np.log10(df1.iloc[:,i]+sf).to_numpy()[:, None]
        pts = np.linspace(-3, 6, 100)[:, None]
        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(P1)
        logkde = kde.score_samples(pts)
        ax.plot(pts,np.exp(logkde))
    # ax.set_ylim([0,0.8])
    ax.legend(df1.columns)
    fig.suptitle('Distribution', fontsize=20)
    plt.xlabel('Log10 TPM', fontsize=14)
    plt.ylabel('Empirical PDF', fontsize=14)



def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df



def stack_fig(df1, position, character, group_names, cols=None, idx_col=None):
    """
    This function takes in panda dataframes with multiple columns and stacks them into one columns with values and adds a columns with group ids...

    Parameters:
        df1 - Pandas DataFrame
            The input dataframe in the form with multiple columns that need to be stacked

        position - int
            An integer that specifies where in the string we want to check what group
        character - string in quotes
            the letter that will differentiate the groups
        group_names - a list in square brackets
            give a list of two names for the group column in the correct order
        
        cols - columns to stack
            if set to None it will use all columns in df

        idx_col - integer
            an integer that will specify which column to reset index. If None it will use first column from df

    """
    
    if cols is None:
        cols = df1.columns

    if idx_col is None:
        idx_col = cols[0]

    df2 = df1.loc[:, cols].set_index(idx_col)
    df3 = df2.T.stack().reset_index().rename(columns={'level_0': 'sample', 0: 'tpm'})
    group_list = [group_names[0] if row['sample'][position].lower() == character else group_names[1] for idx, row in
                  df3.iterrows()]
    df3['group'] = group_list
    return df3








