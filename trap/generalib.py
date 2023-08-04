
import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt


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


    df2=df1.loc[:,cols].set_index(idx_col)
    df3=df2.T.stack().reset_index().rename(columns={'level_0': 'sample', 0: 'tpm'})

    grouplist = [group_names[0] if row['sample'][position] == character else group_names[1] for idx,row in df3.iterrows()]
    df3['group'] = grouplist

    return df3


def kdensity(df1, sf=0.01):

    fig,ax=plt.subplots(figsize=(8, 6))
    for i in range(0,len(df1.columns)):
        P1 = np.log10(df1.iloc[:,i]+sf).to_numpy()[:, None]
        pts = np.linspace(-2, 6, 100)[:,None]
        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(P1)
        logkde = kde.score_samples(pts)
        ax.plot(pts,np.exp(logkde))
    # ax.set_ylim([0,0.8])
    ax.legend(df1.columns)
    fig.suptitle('distribution', fontsize=20)
    plt.xlabel('Log10 TPM', fontsize=14)
    plt.ylabel('Empirical PDF', fontsize=14)


def stack_fig1(df1, position, character1, character2, character3, group_names, cols=None, idx_col=None):
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


    df2=df1.loc[:,cols].set_index(idx_col)
    df3=df2.T.stack().reset_index().rename(columns={'level_0': 'sample', 0: 'tpm'})


    
    grouplist=[]
    for idx,row in df3.iterrows():
        if row['sample'][position] == character1:
            groupname = group_names[0]
        elif row['sample'][position] == character2:
            groupname = group_names[1]
        elif row['sample'][position] == character3:
            groupname = group_names[2]
    
        grouplist.append(groupname)

    df3['group'] = grouplist
    return df3