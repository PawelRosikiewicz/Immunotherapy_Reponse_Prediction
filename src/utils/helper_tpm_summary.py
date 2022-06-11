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



# Function, ......................................
def tpm_summary(df, n=None, plot=False, deg=2, hue=None):
    ''' creates fast summary for gene expression file, with tmp, 
        it helps me to set, or test applied thresholds, for filtering
        
        parameters
        . df; pandas Data Frame
        . n; number of randomly selected columns for creating plots, and summary table
        . plot; bool, if True, it shows the plot (see below)
        . deg; degree for polyft function, used for making trneline on mea_tpm~sd.tpm plot
        
        returns
        . pandas series wiht tmp statistics - stats descibed directly in pd.series, 
        . plot; if plot=True, with 2 subplots, showing % of genes/sample and sample/gene with tmp>0
    
        future developments
        . trendline: https://stackoverflow.com/questions/66040288/python-smoothing-2d-plot-trend-line
            here with plotly: https://plotly.com/python/linear-fits/
    '''
    # ensure, all tmp's are numeric
    arr = df.values.astype("float64")
    
    # hue vector
    if(hue is None):
        hue = [0]*arr.shape[0]
    else:
        pass
    hue_unique = pd.Series(hue).unique()
    
    # subset df columns, eg sample 1000 columns out of 10000
    if(n is None):
        arr = arr
        n = arr.shape[1]
    else:    
        idx = np.random.choice(list(range(df.shape[1])), size=n, replace=False)
        arr = arr[:,idx].copy()
        
    # create summary for each gene/column
    res=[]
    
    # .. one columns for each hue value
    for hue_value in hue_unique:
        arr_sub = arr[hue==hue_value,:] 
        res.append({
            "sample nr":arr_sub.shape[0],
            "gene nr":arr_sub.shape[1],
            "- VALUES -":"",
            "min tpm value recorded": np.round(arr_sub.min(),1),
            "median tpm value recorded": np.round(np.median(arr_sub),1),
            "mean tpm value recorded": np.round(arr_sub.mean(),1),
            "max tpm value recorded": np.round(arr_sub.max(),1),
            "- DETECTED IN -":"",
            "mean % of genes expressed per sample": f'{100-np.round(((arr_sub==0).sum(axis=1)/n*100).mean(),1)}%',
            "% of genes not expressed in any sample":f'{np.round(((arr_sub==0).sum(axis=0)==arr_sub.shape[0]).sum()/n*100,1)}%', # all are zero
            "% of genes expressed in at least 50%  of samples":f'{np.round(((arr_sub==0).sum(axis=0)<=arr_sub.shape[0]/2).sum()/n*100,1)}%', 
          # expressed in at least half of the samples
            "% of genes expressed in all samples":f'{np.round(((arr_sub==0).sum(axis=0)<=0).sum()/n*100,1)}%', # expressed in all samples
        })

    res_df = pd.DataFrame(res, columns=hue_unique) 
    return pd.DataFrame(res)



# Function, ......................................
def tpm_plots(df, n=None, deg=2, hue=None,title=None, color=None):
    ''' creates fast summary for gene expression file, with tmp, 
        it helps me to set, or test applied thresholds, for filtering
        
        parameters
        . df; pandas Data Frame
        . n; number of randomly selected columns for creating plots, and summary table
        . plot; bool, if True, it shows the plot (see below)
        . deg; degree for polyft function, used for making trneline on mea_tpm~sd.tpm plot
        
        returns
        . pandas series wiht tmp statistics - stats descibed directly in pd.series, 
        . plot; if plot=True, with 2 subplots, showing % of genes/sample and sample/gene with tmp>0
    
        future developments
        . trendline: https://stackoverflow.com/questions/66040288/python-smoothing-2d-plot-trend-line
            here with plotly: https://plotly.com/python/linear-fits/
    '''
    # ensure, all tmp's are numeric
    arr = df.values.astype("float64")
    
    # set colors
    if(color is None):
        colors = ["navy", "darkorange", "yellowgreen", "cyan", "gold"]
    else:
        colors = [color]*100
    
    # hue vector
    if(hue is None):
        hue = [0]*arr.shape[0]
    else:
        pass
    hue_unique = pd.Series(hue).unique()
    
    # subset df columns, eg sample 1000 columns out of 10000
    if(n is None):
        arr = arr
        n = arr.shape[1]
    else:    
        idx = np.random.choice(list(range(df.shape[1])), size=n, replace=False)
        arr = arr[:,idx].copy()
        
    # create summary for each gene/column
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(9,3), facecolor="white")
    if(title is not None):
        fig.suptitle(title)
    else:
        pass
    
    # create subplots on each dataset provided with different hue value
    for color, hue_value in zip(colors, hue_unique):
        arr_sub = arr[hue==hue_value,:] 

        # subplot 1. in how many samples, tmp is>0
        nonzero_tmp_per_gene = (1-(arr_sub==0).sum(axis=0)/arr.shape[0])*100
        axs[0].hist(nonzero_tmp_per_gene, alpha=0.5, color=color, label=hue_value, align='mid')
        axs[0].set_xlim(0,100)
        
        # subplot 2
        nonzero_tmp_per_sample = (1-(arr_sub==0).sum(axis=1)/n)*100
        axs[1].hist(nonzero_tmp_per_sample, alpha=0.5, color=color, label=hue_value)
        
        # .. data
        xdata = arr_sub.mean(axis=0)
        ydata = arr_sub.std(axis=0)
        axs[2].scatter(y=ydata, x=xdata, marker='h', alpha=0.5, color=color, label=hue_value)
        
        # ... trendline
        z = np.polyfit(xdata, ydata, deg)
        x_new = np.linspace(0, xdata.max(), 300)
        p = np.poly1d(z)
        axs[2].plot(x_new, p(x_new), color="black")        
        
    # aestetics
    axs[0].set_title("% of samples with detected\nexpression/transcripts")
    axs[0].set_ylabel("frequency")
    axs[0].set_xlabel("% of samples")
    axs[0].grid(lw=0.2)

    # subplot 2. 
    axs[1].set_title("% og genes expressed\n per sample")
    axs[1].set_ylabel("frequency")
    axs[1].set_xlabel("% of genes/sample")
    axs[1].grid(lw=0.2)

    # subplot 3. mean tpm ~ sd
    axs[2].set_title("standard deviation~the mean\nfor each gene")
    axs[2].set_ylabel("sd")
    axs[2].set_xlabel("mean")
    axs[2].grid(lw=0.2)
    
    # add legend
    if len(hue_unique)>1:
        axs[0].legend()
        axs[0].legend()
        axs[0].legend()
    else:
        pass

    sns.despine()
    fig.tight_layout()
    plt.show();