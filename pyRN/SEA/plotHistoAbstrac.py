#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 10:03:54 2022

@author: pmaldona
"""
import numpy as np

def plotHistAbstrac(subplot,abstList,SpIdArray):
    '''
    

    Parameters
    ----------
    subplot : Matplotlib subplot object.
    abstList : List of abstractions.
    SpIdArray : Array of species names.

    Returns
    -------
    Generates a subplot of frequencies of abstractions.

    '''    
    
    SpAbstList=list(map(lambda x: '{'+",".join(SpIdArray[np.where(x)[0]])+'}',abstList))
    
    subplot.hist(SpAbstList)
    subplot.tick_params(labelrotation=45)
    subplot.set_title("Abstraction distribution")
    subplot.set_ylabel("Frecuencies")
    
def abstractionsFromRwSimpleDict(SimpleRWDict, w=None, abstraction_type='c'):
    '''
    In:  SimpleRWDict: Dictonary of simple random walk 
         w: Kyes of SimpleRWDict Dictonary considered, if w=None, it use all the 
         Dictionary
         abstraction_type: corrspond to the type of abstraction to be considered, if 
         its 'all' used both perturbation and convegence states
    Out: list of abstracion of a simple Random Walk
    '''  

    if w is None:
        w=SimpleRWDict.keys()
        
    for i in w:
        
        if abstraction_type=='all':
            
            abs_1 = SimpleRWDict[i]['p']# Converts nested dictionary to 2-dimensional array
            abs_2 = SimpleRWDict[i]['c']# Converts nested dictionary to 2-dimensional array 
            abstractions = [item for sublist in zip(abs_1, abs_2) for item in sublist]
        else:
            abstractions = SimpleRWDict[i][abstraction_type]# Converts nested dictionary to 2-dimensional array
            
        return abstractions
    
    
    