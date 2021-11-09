#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from CRNS import CRNS

#loading form sbml file 
file="./networks/BIOMD0000000013.xml"
RN = CRNS.from_sbml(file)

#loading form a text file please refer to rn_test.txt to see example
file="./networks/rn_test.txt"
RN = CRNS.form_txt(file)

print(RN.mp)
print(RN.mr)
print(RN.sp)
print(RN.sp_n)

# Considering an small set
sp_set = RN.sp[[5,7,9]]
print(sp_set)

# geretaed closure for a given set
close_sp_set=RN.closure(sp_set)
print(close_sp_set)

# One step generated closure
one_clos_sp = RN.clos_one(sp_set)
print(one_clos_sp)

# reactive semi-self-mantaining
print(RN.is_ssm(sp_set))
print(RN.is_ssm(close_sp_set))

# stochomerically semi-self-mantaining
print(RN.is_s_ssm(sp_set))
print(RN.is_s_ssm(close_sp_set))

# if subset is self-mantainig
print(RN.is_sm(sp_set))
print(RN.is_sm(close_sp_set))

# Generation of reaction equivalence clases and basic sets
RN.gen_basics()
 
