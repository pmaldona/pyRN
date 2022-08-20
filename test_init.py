#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:34:08 2022

@author: pmaldona
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from pyRN import pyRN

# loading form a text file please refer to rn_test.txt to see example,
# textfiles preserve antimony general strucutre in terms of reaction and 
# reaction arrows.
file="networks/rn_test.txt"
RN = pyRN.from_txt(file)

# Alternative sbml files can be loaded. 
# file="networks/PW000035.sbml"
# RN = pyRN.from_sbml(file,False)

# Basic Variables can easily obtain
print("Vector of species")
print(RN.sp)
print("Vector of species names, in case of an initialization form an smbl file result can be different")
print(RN.sp_n)
print("Display species")
print(RN.sin_print_sp())
print("Display reactions")
print(RN.sin_print_r())
print("Reactive stochiometric matrix")
print(RN.mr)
print("Productive stoichimetric matrix")
print(RN.mp)
print("Stoichiometric Matrix")
print(RN.mp-RN.mr)
print("Inflow species")
print(RN.is_inflow(RN.sp,True))
print("Outflow species")
print(RN.is_outflow(RN.sp,True))
RN.plot_S()
nt=RN.display_RN()
nt.show("RN.html")