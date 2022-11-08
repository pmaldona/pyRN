#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:53:50 2022

@author: pmaldona
"""


from .RNSRW import RNSRW
from .SAE.plot_raw import plot_raw
from .SAE.plot_change import plot_change
from .SAE.plot_complexitychange import plot_complexitychange
from .SAE.plot_hasse import plot_hasse
from .SAE.plot_hasse import plot_hasse_convergence_and_perturbation
from .SAE import sos 

class RNSEA(RNSRW):
    
    # Raw plot inherence for SAE scripts (Simon development). Recive as input 
    # a list of abstractions, subplot matplotlib pyplot element and a string
    # as title
    def plotRawRw(self,subplot, abstractions, title=''):
        plot_raw(subplot, abstractions, title=title)
        
    # Raw plot inherence for SAE scripts (Simon development). Recive as input 
    # Walk index of the SimpleRw Dictionary, and the type of abstraction ('c' 
    # for convergent states and 'p' for perturbated stares), title and subplot same as function before
    def plotRawSimpleRw(self,subplot, walk_index=0, abstractions_type='c', title=''):
        plot_raw(subplot, self.RwSimpleListDictDf[walk_index][abstractions_type], title=title)
    
    # Function that displays the evolution of an Random Walk in terms of growment and contractions
    def plotChangeRw(self,subplot, abstractions, change_function=sos.normalized_hamming_distance, title='', show_indices=True, index_spacing=1, legend=True):
        '''
        In:  mandatory:
             subplot (matplotlib.pyplot subplot)
             abstractions (list), list of binary arrays that represent the abstractions in the random walk

             optional:
             change_function (function), function that calculates the distance between pairs of abstractions
             show_indices (boolean), chooses if indices of transitions will be displayed
             index_spacing (float), modifies spacing between indices
             legend (boolean), chooses if a legend will be displayed
        '''
        
        plot_change(subplot, abstractions, change_function=change_function, title=title, 
                    show_indices=show_indices, index_spacing=index_spacing, legend=legend)
    # Variation of latter funtion but for a simple random walk. 
    def plotChangeSimpleRw(self,subplot, walk_index=0, abstractions_type='c', change_function=sos.normalized_hamming_distance, title='', show_indices=True, index_spacing=1, legend=True):
        '''
        In:  mandatory:
             subplot (matplotlib.pyplot subplot)
             walk_index, index of the random walk to be taken into account
             abstraction_type, type of abstraction 'c' for convergent and 'p' for perturbated 

             optional:
             change_function (function), function that calculates the distance between pairs of abstractions
             show_indices (boolean), chooses if indices of transitions will be displayed
             index_spacing (float), modifies spacing between indices
             legend (boolean), chooses if a legend will be displayed
        '''
        
        plot_change(subplot, self.RwSimpleListDictDf[walk_index][abstractions_type], change_function=change_function, title=title, 
                    show_indices=show_indices, index_spacing=index_spacing, legend=legend)
    
    # Plot funtion of the random walk in the hasse diagram.
    def plotHasseRw(self, subplot, abstractions, title='', loga=False):
        '''
        In:  subplot (matplotlib.pyplot subplot)
             abstractions (list), list of binary arrays that represent the abstractions in the random walk
             change_function (function), function that calculates the distance between pairs of abstractions
            loga (boolean), logarithmic scale on x_axis (recommended for larger number of species
        '''
        plot_hasse(subplot, abstractions, title=title, loga=loga)
        
    # Variant of the latter function but considering the simple random walks.
    def plotHasseSimpleRw(self, subplot, walk_indexes=None, abstractions_type='c', title='', loga=False):
        '''
        In:  subplot (matplotlib.pyplot subplot)
             walk_indexes, list of indexes of the random walks to be taken into account
             abstraction_type, type of abstraction 'c' for convergent and 'p' for perturbated 
             change_function (function), function that calculates the distance between pairs of abstractions
             loga (boolean), logarithmic scale on x_axis (recommended for larger number of species
        '''
        if walk_indexes is None:
            walk_indexes=self.RwSimpleListDictDf.keys()
            
        abstractions=[]
        for i in walk_indexes:
            abstractions.append(self.RwSimpleListDictDf[i][abstractions_type])
        plot_hasse(subplot, abstractions, title=title, loga=loga)
        
    # Plot funtion of the random walk in the hasse diagram considering convergent and perturbed states   
    def plotHasseConvergenceAndPerturbationRw(self, subplot, convergent_abstractions, perturbation_abstractions, title='', loga=False):
        '''
        In:  subplot (matplotlib.pyplot subplot)
             abstractions (list), list of binary arrays that represent the abstractions in the random walk
             change_function (function), function that calculates the distance between pairs of abstractions
            loga (boolean), logarithmic scale on x_axis (recommended for larger number of species
        '''
        plot_hasse_convergence_and_perturbation(subplot, convergent_abstractions, perturbation_abstractions, title=title, loga=loga)
    
    # Variant of the latter function but considering a simple random walk  
    def plotHasseConvergenceAndPerturbationSimpleRw(self, subplot, walk_index=0, title='', loga=False):
        '''
        In:  subplot (matplotlib.pyplot subplot)
             walk_index, index of the random walk to be taken into account
             change_function (function), function that calculates the distance between pairs of abstractions
            loga (boolean), logarithmic scale on x_axis (recommended for larger number of species
        '''
        
        plot_hasse_convergence_and_perturbation(subplot, self.RwSimpleListDictDf[walk_index]['c'], 
                                                self.RwSimpleListDictDf[walk_index]['p'], title=title, loga=loga)
    def plotComplexityChangeRw(subplot, abstractions, abstractions_df, title='', show_indices=True, index_spacing=1):
        '''
        In:  mandatory:
             subplot (matplotlib.pyplot subplot)
             abstractions (list), list of binary arrays that represent the abstractions in the random walk
             abstractions_df (pandas dataframe)

             optional:
             change_function (function), function that calculates the distance between pairs of abstractions
             show_indices (boolean), chooses if indices of transitions will be displayed
             index_spacing (float), modifies spacing between indices
             legend (boolean), chooses if a legend will be displayed

        '''
        plot_complexitychange(subplot, abstractions, abstractions_df, title=title, 
                              show_indices=show_indices, index_spacing=index_spacing)
    
    def plotComplexityChangeSimpleRw(subplot, abstractions, abstractions_df, title='', show_indices=True, index_spacing=1):
        '''
        In:  mandatory:
             subplot (matplotlib.pyplot subplot)
             abstractions (list), list of binary arrays that represent the abstractions in the random walk
             abstractions_df (pandas dataframe)

             optional:
             change_function (function), function that calculates the distance between pairs of abstractions
             show_indices (boolean), chooses if indices of transitions will be displayed
             index_spacing (float), modifies spacing between indices
             legend (boolean), chooses if a legend will be displayed

        '''
        plot_complexitychange(subplot, abstractions, abstractions_df, title=title, 
                              show_indices=show_indices, index_spacing=index_spacing)    
   