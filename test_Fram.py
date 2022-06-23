#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 18:19:05 2022

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from CRNS import CRNS
from bitarray import frozenbitarray
import networkx as nx
from bitarray import bitarray as bt


#loading form sbml file 
# file="networks/BIOMD0000000013.xml"
file="./networks/Farm_sbml.xml"
RN = CRNS.from_sbml(file,False)

# loading form a text file please refer to rn_test.txt to see example
# file="networks/rn_test.txt"
# RN = CRNS.form_txt(file)

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
print(RN.is_ssm(RN.sp))

# stochomerically semi-self-mantaining
print(RN.is_s_ssm(sp_set))
print(RN.is_s_ssm(close_sp_set))
print(RN.is_s_ssm(RN.sp))

# if subset is self-mantainig
print(RN.is_sm(sp_set))
print(RN.is_sm(close_sp_set))
print(RN.is_sm(RN.sp))


# Generation of reaction equivalence clases and basic sets
RN.gen_basics()
print(RN.sp_b)
print(RN.conn)
print(RN.dyn_conn)

# Generation of the synergic structure
RN.gen_syn_str()
print(len(RN.syn_str.nodes()))
print(len(RN.syn_str.edges()))
print(len(RN.ssms))
ssm_nodes = [x for x,y in RN.syn_str.nodes(data=True) if y['is_ssm']]

# Generation of the semi-self-maintained structure
RN.gen_ssm_str()
print(len(RN.ssm_str.nodes()))
print(len(RN.ssm_str.edges()))
print(len(RN.ssms))


