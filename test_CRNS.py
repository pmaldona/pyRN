#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from CRNS import CRNS
from bitarray import frozenbitarray
import networkx as nx

#loading form sbml file 
# file="networks/BIOMD0000000013.xml"
file="networks/PW000035.sbml"
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
print(RN.syn_str.nodes())
print(RN.syn_str.edges())
print(RN.ssms)

selected_edges = [(u,v,e) for u,v,e in RN.syn_str.edges(data=True) if e['syn'] == True]

print(RN.syn_str.nodes[frozenbitarray('000101')])
for u,v,e in RN.syn_str.edges(data=True):
    if e['syn']==True:
        print ([u,v,RN.syn_str[u][v]])

# Generation of the synergic structure
RN.gen_ssm_str()
print(RN.ssm_str.nodes())
print(RN.ssm_str.edges())
print(RN.ssms)

# generation of minimal generators
RN.gen_mgen()
print(RN.mgen)

# Generation of all proto synergies
RN.all_syn()   
print(RN.syn)
print(RN.syn_p)
