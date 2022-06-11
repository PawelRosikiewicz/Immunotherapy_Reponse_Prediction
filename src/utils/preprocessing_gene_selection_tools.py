# ********************************************************************************** #
#                                                                                    #
#   Project: Ardigen task                                                            #                        
#   Author: Pawel Rosikiewicz                                                        #
#   Contact: prosikiewicz(a)gmail.com                                                #
#                                                                                    #
#   License: MIT License                                                             #
#   Copyright (C) 2022.06.04 Pawel Rosikiewicz                                       #
#                                                                                    #
# Permission is hereby granted, free of charge, to any person obtaining a copy       #
# of this software and associated documentation files (the "Software"), to deal      #
# in the Software without restriction, including without limitation the rights       #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell          #
# copies of the Software, and to permit persons to whom the Software is              #
# furnished to do so, subject to the following conditions:                           #
#                                                                                    # 
# The above copyright notice and this permission notice shall be included in all     #
# copies or substantial portions of the Software.                                    #
#                                                                                    #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR         #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,           #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE        #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER             #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,      #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE      #
# SOFTWARE.                                                                          #
#                                                                                    #
# ********************************************************************************** #



#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os # allow changing, and navigating files and folders, 
import sys
import re # module to use regular expressions, 
import glob # lists names in folders that match Unix shell patterns
import random # functions that use and generate random numbers

import numpy as np # support for multi-dimensional arrays and matrices
import pandas as pd # library for data manipulation and analysis
import seaborn as sns # advance plots, for statistics, 
import matplotlib as mpl # to get basic plt   functions, heping with plot mnaking 
import matplotlib.pyplot as plt # for making plots, 
import scipy.stats as stats




# Function, ..............................................
def compare_gene_expression(x, y, method="median"):
    ''' simple funciton that calulates gene expression foldchnage between two classes of samples
        and ttstudent test, then it data frame, with list of all genes, sorted with pvalue from ttest, 
        from the most, to the least signifficant,
        . x, y - dataframe and pd. series repsectively, with the data, and labels, 
        . method - str, "median", of "mean"
        ...............
        important: ttest shoudl be used only with relatively large number of replicates in each gorup, 
        eg >20, For more informaiton see these articles:
           https://academic.oup.com/ije/article/39/6/1597/736515
           https://doi.org/10.1016/j.csbj.2021.05.040
    '''
    # test input df, & work on df copy,
    assert type(x) == pd.DataFrame, "x_train Incorrect obj type"
    assert type(y) == pd.Series, "y_train Incorrect obj type"
    x, y = x.copy(), y.copy()
        
    # divide the set into two group
    ttest=[]
    Log2FC=[]
    for idx in range(x.shape[1]):
        one_row = x_transf.iloc[:,idx].values
        a = one_row[y==0]
        b = one_row[y==1]

        # .. ttest
        ttest.append((stats.ttest_ind(a, b).pvalue))

        # Log2FC
        if method=="median": Log2FC.append(np.median(a)/np.median(b)-1)
        if method=="mean": Log2FC.append(np.mean(a)/np.mean(b)-1)
            
    # store results in nice dataframe
    results = pd.DataFrame([ttest,Log2FC]).transpose()
    results.columns = ['Pval', 'Log2FC']
    results.index = x.columns   
    
    # calculate LogPval
    results['LogPval'] = -(np.log(results.loc[:,"Pval"]))
    
    return results.sort_values(by="Pval", ascending=True)

  
  
# Function, ..............................................
def select_genes_and_make_volcano_plot(
    df, xname="Log2FC", yname="LogPval", pname="Pval", 
    title="Volcano Plot", figsize=(10,4),Ptr=0.05, Log2FCtr=2, create_plot=True
):
    ''' This fucntions, allows selecting differencially expressed 
        genes with PValue, and Log2 fold chnage thresholds, 
        it returns, the table wiht selected genes, (gene names are indexes),
        and if create_plot is True, it also returns volcano polot wiht selected genes, 
        and basic informaiton on their number
    
        parameters:
        . df; Pandas datagrame, input data
        . xname; def. "Log2FC"; column name wiht log 2 fold chnage data in df, also y-axis on a plot
        . yname; def."LogPval"; column name wiht -log10(p-value) data in df, also y-axis on a plot
        . pname; def."Pval"; column name wiht p-value data in df, used for adding colors, 
            on a scatter points, and selection fo the points for the table
        . title; str, added on top of the title on a volcano plot,
        . figsize, tuple, with two int,
        . Ptr; float, threshold on p-value 
        . Log2FCtr; float or int, threshold on Log2FC 
        . create_plot; if True, returns volcano plot       

    '''
    # test input df, & work on df copy,
    assert type(df) == pd.DataFrame, "x_train Incorrect obj type"
    df = df.copy()   
    

    # (A) prepare the data

    #.select upregulated & downregulatzed genes
    downregulated_genes = df.loc[(df.loc[:,xname]<(-Log2FCtr)) & (df.loc[:,pname]<=Ptr), :]
    upregulated_genes = df.loc[(df.loc[:,xname]>Log2FCtr) & (df.loc[:,pname]<=Ptr), :]    
    
    
    # (B) volcano plot
    
    if create_plot==True:

        # make figure
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize)

        # scatter plot - all points, 
        ax.scatter( x=df.Log2FC, y=df.LogPval, 
                    marker='o', alpha=1, s=6, color="grey")

        # add  subpltitle  with info on pre-selected postions
        cut_offs = f'cut-off of Pval={Ptr}, and Log2FC=±{Log2FCtr}' 
        position_nr = f'DOWNREGULATED: {downregulated_genes.shape[0]}, UPREGULATED:{upregulated_genes.shape[0]}'
        fig.suptitle(f'{title}\n{cut_offs}\n{position_nr}')       

        # upregulated and downregulated points wiht red and blue colors
        ax.scatter(
            x=downregulated_genes.loc[:,xname], 
            y=downregulated_genes.loc[:,yname], 
                    marker='h', alpha=1, s=5, color="blue")    
        ax.scatter(
            x=upregulated_genes.loc[:,xname], 
            y=upregulated_genes.loc[:,yname], 
                    marker='h', alpha=1, s=5, color="red")   

        # set limits, to make the plot more easy to read
        xlimits = np.abs([downregulated_genes.loc[:,xname].min(),upregulated_genes.loc[:,xname].max()])
        ax.set_xlim(-xlimits.max(), xlimits.max())    

        # horizontal and vertical lines with Pval and Log2FC thresholds
        ax.axhline(-(np.log(Ptr)), lw=0.5, ls="--", color="black")
        ax.axvline(-(Log2FCtr), lw=0.5, ls="--", color="black")
        ax.axvline(Log2FCtr, lw=0.5, ls="--", color="black")

        # ax aestetics
        ax.grid(lw=0.3)
        ax.set_xlabel("Log2FC")
        ax.set_ylabel("-Log10(Pval)")

        # fig layout
        fig.subplots_adjust(top=0.7)
        sns.despine()    
        plt.show();
    
    else:
        pass
    
    
    # (C) return table
    selected_genes = pd.concat([downregulated_genes, upregulated_genes], axis=0)
    return selected_genes