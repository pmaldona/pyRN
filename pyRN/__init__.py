#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 19:15:03 2022

@author: pmaldona

General Reaction Network Class
"""

from .RNLI import RNLI
from .RNSEA import RNSEA
import copy
    
class pyRN(RNLI,RNSEA):
    
    def copy(self):
        return pyRN(copy.copy(self))
    
    pass


    
    
        

