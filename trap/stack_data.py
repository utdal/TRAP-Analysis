#functions used in TRAP analysis

import pandas as pd
import numpy as np
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


tempdf3= pd.read_csv(r'path/to/file')
new_df= stack_fig(tempdf3, position=4, character="N", group_names=['Input', 'IP'], cols=None, idx_col=None)
