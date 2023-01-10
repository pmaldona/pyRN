#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 19:15:03 2022

@author: pmaldona

General Reaction Network Class
"""

from .RNLI import RNLI
from .RNSEA import RNSEA
import pickle
import copy
import inspect
    
class pyRN(RNLI,RNSEA):
    
    def copy(self):
        '''
        

        Returns
        -------
        pyRN objet
            hard copy of pyRN object.

        '''
        return pyRN(copy.copy(self))
    
    pass
        
    def setFromPkl(self,file):
        '''
        

        Parameters
        ----------
        file : pkl. 
            Pikle file contiaing pyRN object.

        Returns
        -------
        None. Copy the variables to the object itself.

        '''
        
        
        # loads the pickle file in a objet
        with open(file, 'rb') as f:
            obj = pickle.load(f)
        
        self.__dict__.update({k: v for k, v in obj.items() if not inspect.isfunction(v)})
        
    
    def saveToPkl(self,file="pyRN_object.pkl"):
        '''
        

        Parameters
        ----------
        file : str, optional
            Pickle file name where the pyRN object will be saved. The default is "pyRN_object.pkl".

        Returns
        -------
        None.

        '''
        # creates a dictionary fo all avriables of the module
        obj = {k: v for k, v in self.__dict__.items() if not k.startswith('__')}

        # Guarda el diccionario en un archivo pickle
        with open(file, 'wb') as f:   
            pickle.dump(obj, f)
            
        
            
        