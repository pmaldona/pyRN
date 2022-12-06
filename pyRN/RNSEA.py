#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 09:53:50 2022

@author: pmaldona
"""

from .RNSRW import RNSRW
from .SEA.plot_raw import plot_raw #./SEA/plot_raw.py import me plot_raw function
from .SEA.plot_change import plot_change
from .SEA.plot_complexitychange import plot_complexitychange
from .SEA.plot_hasse import plot_hasse
from .SEA.plot_hasse import plot_hasse_convergence_and_perturbation
from .SEA.dataframes import dataframesFromLists
from .SEA.plot_markov import plot_markov
from .SEA.plotHistoAbstrac import plotHistAbstrac
from .SEA import sos 
import itertools


class RNSEA(RNSRW):
    

    def plotRawRw(self,subplot, walk_type='simple', walk_index=0, abstraction_type='c', title=''):
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_index : int, optional
            walk of the radom walk repetition to be plot. The default is 0.
        abstraction_type : string, optional
            Typr abstraction selected to be polt. The default is 'c'.
        title : string, optional
            Title of the plot.  The default is ''.

        Returns
        -------
        Changes de subplot displaying the raw data of the random walk.

        '''
        
        plot_raw(subplot, self.RwDict[walk_type][walk_index][abstraction_type].transpose(), title=title)
    

    def plotChangeRw(self, subplot,  walk_type='simple', walk_index=0, abstraction_type='c', change_function=sos.normalized_hamming_distance, title='', show_indices=True, index_spacing=1, legend=True):
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_index : int, optional
            Index of the random walk to be taken into account. The default is 0.
        abstraction_type : string, optional
            Type of abstraction. The default is 'c'.
        change_function : function, optional
            Function that calculates the distance between pairs of abstractions. The default is sos.normalized_hamming_distance.
        title : TYPE, optional
            DESCRIPTION. The default is ''.
        show_indices : TYPE, optional
            DESCRIPTION. The default is True.
        index_spacing : TYPE, optional
            Modifies spacing between indices. The default is 1.
        legend : bool, optional
            Chooses if a legend will be displayed. The default is True.

        Returns
        -------
        Changes de subplot displaying the evolution of an Random Walk in terms of growment and contractions
        '''
        
        plot_change(subplot, self.RwDict[walk_type][walk_index][abstraction_type].transpose().values.tolist(), 
                    change_function=change_function, title=title, 
                    show_indices=show_indices, index_spacing=index_spacing, legend=legend)


    # Variant of the latter function but considering the simple random walks.
    def plotHasseRw(self, subplot, walk_type='simple', walk_indexes=None, abstraction_type='c', title='', loga=False):
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_indexs : range, optional
            choosen walks of the random walk repetition to be plot. The default is None choosing all.
        abstraction_type : string, optional
            Typr abstraction selected to be polt. The default is 'c'.
        title : string, optional
            Title of the plot.  The default is ''.
        loga : bool, optional
            Logarithmic scale on x_axis (recommended for larger number of species). The default is False.

        Returns
        -------
        changes subplot displaying a plot of the random walk in the hasse diagram.

        '''
        
        if walk_indexes is None:
            walk_indexes=self.RwDict[walk_type].keys()
            
        abstractions=[]
        for i in walk_indexes:
            abstractions.append(self.RwDict[walk_type][i][abstraction_type].transpose().values.tolist())
        plot_hasse(subplot, abstractions, title=title, loga=loga)
            
    
    def plotHasseConvergenceAndPerturbationRw(self, subplot, walk_type='simple', walk_index=0,  
                                              convergent_abstraction_type='c',
                                              title='', loga=False):
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_index : int, optional
            walk of the radom walk repetition to be plot. The default is 0.
        convergent_abstraction_type : TYPE, optional
            Convergence abstraction after the simulation. The default is 'c'.
        title : string, optional
            Title of the plot.  The default is ''.
        loga : bool, optional
            Logarithmic scale on x_axis (recommended for larger number of species). The default is False.

        Returns
        -------
        Plot funtion of the random walk in the hasse diagram considering convergent and perturbed states.

        '''
        if walk_type=='mak':
            perturbation_abstractions=(self.RwDict[walk_type][walk_index]['p']>0).transpose().values.tolist()
        else:
            perturbation_abstractions=self.RwDict[walk_type][walk_index]['p'].transpose().values.tolist()
        
        plot_hasse_convergence_and_perturbation(subplot, 
                                                self.RwDict[walk_type][walk_index][convergent_abstraction_type].transpose().values.tolist(), 
                                                perturbation_abstractions, 
                                                title=title, loga=loga)
    
    def getAbstrationTransitionDf(self,walk_type='simple',walk_indexes=None,abstraction_type='c',complexity_type='cc'):
        '''
        

        Parameters
        ----------
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_indexs : range, optional
            choosen walks of the random walk repetition to be plot. The default is None choosing all.
        abstraction_type : string, optional
            Type abstraction selected to be plot. The default is 'c'.
        complexity_type : string, optional
            Type complexity selected to be plot. The default is 'cc'.


        Returns
        -------
        abstractions_df : Pandas DataFrame
            Abstraction dataframe whit resilient data.
        transitions_df : Pandas DataFrame
            Transition dataframe whit resilient data.
        '''
        
        if walk_indexes is None:
            walk_indexes=self.RwDict[walk_type].keys()
            
        abstractions=[]
        complexities=[]
        for i in walk_indexes:
            abstractions.append(self.RwDict[walk_type][i][abstraction_type].transpose().astype(int).values.tolist())
            complexities.append(self.RwDict[walk_type][i][complexity_type])
            
        return dataframesFromLists(abstractions,complexities)
    
    
    def plotMarkovRw(self, subplot, walk_type='simple', walk_indexes=None, abstraction_type='c', complexity_type='cc'):
        
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_indexs : range, optional
            choosen walks of the random walk repetition to be plot. The default is None choosing all.
        abstraction_type : string, optional
            Type abstraction selected to be plot. The default is 'c'.
        complexity_type : string, optional
            Type complexity selected to be plot. The default is 'cc'.


        Returns
        -------
        Plot funtion that modefies subplot displaying the Markovian probabilities of a random walk
        '''
        abstractions_df, transitions_df = self.getAbstrationTransitionDf(walk_type=walk_type, 
                                                                walk_indexes=walk_indexes, 
                                                                abstraction_type=abstraction_type,
                                                                complexity_type=complexity_type)
        plot_markov(subplot,abstractions_df,transitions_df)
    
    def plotHistAbstRw(self,subplot,walk_type='simple',walk_indexes=None, abstraction_type='c'):
        '''
        

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Typo of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_indexs : range, optional
            choosen walks of the random walk repetition to be plot. The default is None choosing all.
        abstraction_type : string, optional
            Type abstraction selected to be plot. The default is 'c'.
        complexity_type : string, optional
            Type complexity selected to be plot. The default is 'cc'.


        Returns
        -------
        changes subplot, plot of histogram of current abstractions.

        '''

        if walk_indexes is None:
            walk_indexes=self.RwDict[walk_type].keys()
            
        abstractions=[]
        for i in walk_indexes:
            abstractions.append(self.RwDict[walk_type][i][abstraction_type].transpose().values.tolist())
        
        abstractions = list(itertools.chain(*abstractions))
        plotHistAbstrac(subplot, abstractions, self.SpIdStrArray)
        
    def plotComplexityChangeRw(self, subplot, walk_type='simple', walk_index=0, 
                               abstraction_type='c', complexity_type='cc', 
                               title='', show_indices=True, index_spacing=1, legend=True):   
        '''
         

        Parameters
        ----------
        subplot : matplotlib subplot object
            Subplot where the plot will be pointed.
        walk_type : string, optional
            Type of randomwak selected form the self.RwDict dictonary. The default is 'simple'.
        walk_index : int, optional
            walk of the radom walk repetition to be plot. The default is 0.
        abstraction_type : string, optional
            Type abstraction selected to be polt. The default is 'c'.
        complexity_type : string, optional
            Type complexity selected to be plot. The default is 'cc'.
        title : string, optional
            Title of the plot. The default is ''.
        show_indices : bool, optional
            If the indexes of the iteration are showed. The default is True.
        index_spacing : TYPE, optional
            Modifies spacing between indices. The default is 1.
        legend : bool, optional
            Chooses if a legend will be displayed. The default is True.

        Returns
        -------
            Plot the growth and contraction in complexity in a evolutornay Simple random walk.
        '''
        

        
        abstractions = self.RwDict[walk_type][walk_index][abstraction_type].transpose().astype(int).values.tolist()   
        abstractions_df, transitions_df = self.getAbstrationTransitionDf(walk_type=walk_type, 
                                                                walk_indexes=[walk_index], 
                                                                abstraction_type=abstraction_type,
                                                                complexity_type=complexity_type)
        
        plot_complexitychange(subplot, abstractions, abstractions_df, title=title, 
                              show_indices=show_indices, index_spacing=index_spacing)
        