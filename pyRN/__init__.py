#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 19:15:03 2022

@author: pmaldona

General Reaction Network Class
"""

from .RNDS import RNDS
from .RNSRW import RNSRW
from .CRNS import CRNS
import copy
    
class pyRN(RNSRW,CRNS,RNDS):
    
    def copy(self):
        return pyRN(copy.copy(self))
    
    pass


    
    
        

